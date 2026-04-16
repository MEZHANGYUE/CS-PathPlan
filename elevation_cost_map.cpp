#include "elevation_cost_map.hpp"
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <iomanip>

#ifdef HAVE_GDAL
#include <cpl_string.h>
#include <vector>
#include <limits>
#endif

ElevationCostMap::ElevationCostMap() {}
ElevationCostMap::~ElevationCostMap() = default;

bool ElevationCostMap::loadElevationData(const std::string &path)
{
  std::string path_to_load = path;
  struct stat stat_buf;
  if (stat(path.c_str(), &stat_buf) == 0) {
      long long size = stat_buf.st_size;
      long long limit = 200LL * 1024 * 1024; // 200MB
      if (size > limit) {
          std::string ovr_path = path + ".ovr";
          struct stat stat_ovr;
          if (stat(ovr_path.c_str(), &stat_ovr) == 0) {
              std::cout << "Elevation file is large (" << size << " bytes). Found .ovr file, using it instead: " << ovr_path << std::endl;
              path_to_load = ovr_path;
          }
      }
  }

#ifdef HAVE_GDAL
  GDALAllRegister();
  std::cerr << "ElevationMap: opening TIFF: " << path_to_load << std::endl;
  GDALDataset *poDS = (GDALDataset*)GDALOpen(path_to_load.c_str(), GA_ReadOnly);
  if (!poDS) {
    std::cerr << "ElevationMap: GDALOpen failed for " << path_to_load << "\n";
    return false;
  }
  int full_w = poDS->GetRasterXSize();
  int full_h = poDS->GetRasterYSize();
  double gt[6];
  if (poDS->GetGeoTransform(gt) == CE_None) {
    elev_geo_.origin_x = gt[0];
    elev_geo_.pixel_w = gt[1];
    elev_geo_.rot_x = gt[2];
    elev_geo_.origin_y = gt[3];
    elev_geo_.rot_y = gt[4];
    elev_geo_.pixel_h = gt[5];
    std::cerr << "ElevationMap: geotransform: origin(" << elev_geo_.origin_x << ", " << elev_geo_.origin_y << ") pixel_w=" << elev_geo_.pixel_w << " pixel_h=" << elev_geo_.pixel_h << std::endl;
  } else {
    elev_geo_.origin_x = 0.0; elev_geo_.pixel_w = 1.0; elev_geo_.rot_x = 0.0;
    elev_geo_.origin_y = 0.0; elev_geo_.rot_y = 0.0; elev_geo_.pixel_h = -1.0;
  }

  const uint64_t bytes_needed = static_cast<uint64_t>(full_w) * static_cast<uint64_t>(full_h) * sizeof(float);
  const uint64_t TARGET_BYTES = 200ULL * 1024ULL * 1024ULL;

  GDALRasterBand *band = poDS->GetRasterBand(1);
  int existing_ov = band ? band->GetOverviewCount() : 0;
  bool found_suitable_ov = false;
  if (existing_ov > 0) {
    for (int i = 0; i < existing_ov; ++i) {
      GDALRasterBand *ob = band->GetOverview(i);
      if (ob) {
        uint64_t ow = static_cast<uint64_t>(ob->GetXSize());
        uint64_t oh = static_cast<uint64_t>(ob->GetYSize());
        if (ow * oh * sizeof(float) <= TARGET_BYTES) {
          found_suitable_ov = true;
          break;
        }
      }
    }
  }

  if (bytes_needed > TARGET_BYTES && !found_suitable_ov) {
     bool ok = performDownsampling(poDS, full_w, full_h, bytes_needed, TARGET_BYTES, path_to_load);
     GDALClose(poDS);
     if (ok) {
       printElevationInfo();
       return true;
     } else {
       return false;
     }
  }

  if (bytes_needed > TARGET_BYTES && existing_ov == 0) {
    std::cerr << "ElevationMap: large raster and no overviews; attempting to build external overviews..." << std::endl;
    GDALClose(poDS);
    CPLSetConfigOption("GDAL_TIFF_OVR", "YES");
    GDALDataset *poDSup = (GDALDataset*)GDALOpen(path_to_load.c_str(), GA_Update);
    if (poDSup) {
      std::vector<int> ov_levels;
      int w = full_w, h = full_h;
      int level = 2;
      while ((w / level) >= 25600 || (h / level) >= 25600) {
        ov_levels.push_back(level);
        level *= 2;
      }
      if (!ov_levels.empty()) {
        int nBands = 1;
        int panBands[1] = {1};
        const char *resampling_try = "MAX";
        CPLErr berr = poDSup->BuildOverviews(resampling_try, static_cast<int>(ov_levels.size()), ov_levels.data(), nBands, panBands, nullptr, nullptr);
        if (berr != CE_None) {
          resampling_try = "AVERAGE";
          berr = poDSup->BuildOverviews(resampling_try, static_cast<int>(ov_levels.size()), ov_levels.data(), nBands, panBands, nullptr, nullptr);
        }
        if (berr != CE_None) {
          std::cerr << "ElevationMap: BuildOverviews failed for " << path_to_load << "\n";
        } else {
          std::cerr << "ElevationMap: external overviews built for " << path_to_load << " using resampling=" << resampling_try << "\n";
        }
      }
      GDALClose(poDSup);
    } else {
      std::cerr << "ElevationMap: cannot open " << path_to_load << " in update mode to build overviews; continuing" << std::endl;
    }
    CPLSetConfigOption("GDAL_TIFF_OVR", nullptr);
    poDS = (GDALDataset*)GDALOpen(path_to_load.c_str(), GA_ReadOnly);
    if (!poDS) {
      std::cerr << "ElevationMap: GDALOpen failed for " << path_to_load << " after overview attempt\n";
      return false;
    }
    band = poDS->GetRasterBand(1);
    existing_ov = band ? band->GetOverviewCount() : 0;
  }

  GDALRasterBand *readBand = band;
  int use_ov_index = -1;
  if (band && existing_ov > 0) {
    for (int i = existing_ov - 1; i >= 0; --i) {
      GDALRasterBand *ob = band->GetOverview(i);
      if (!ob) continue;
      uint64_t ow = static_cast<uint64_t>(ob->GetXSize());
      uint64_t oh = static_cast<uint64_t>(ob->GetYSize());
      uint64_t need = ow * oh * sizeof(float);
      if (need <= TARGET_BYTES) {
        readBand = ob;
        use_ov_index = i;
        elev_width_ = static_cast<int>(ow);
        elev_height_ = static_cast<int>(oh);
        double scale_x = static_cast<double>(full_w) / static_cast<double>(ow);
        double scale_y = static_cast<double>(full_h) / static_cast<double>(oh);
        elev_geo_.pixel_w = elev_geo_.pixel_w * scale_x;
        elev_geo_.pixel_h = elev_geo_.pixel_h * scale_y;
        break;
      }
    }
  }

  if (use_ov_index == -1) {
    elev_width_ = full_w;
    elev_height_ = full_h;
  }

  elev_data_.assign(static_cast<size_t>(elev_width_) * static_cast<size_t>(elev_height_), std::numeric_limits<float>::quiet_NaN());
  CPLErr err;
  if (readBand) {
    err = readBand->RasterIO(GF_Read, 0, 0, elev_width_, elev_height_, elev_data_.data(), elev_width_, elev_height_, GDT_Float32, 0, 0);
  } else {
    err = CE_Failure;
  }
  if (err != CE_None) {
    std::cerr << "ElevationMap: GDAL RasterIO failed\n";
    GDALClose(poDS);
    return false;
  }
  GDALClose(poDS);
  elev_valid_ = true;
  std::cerr << "ElevationMap: finished reading raster (width=" << elev_width_ << " height=" << elev_height_ << ")" << std::endl;
  return true;
#else
  std::cerr << "ElevationMap: GDAL not available at build time. Cannot load TIFF." << std::endl;
  return false;
#endif
}

bool ElevationCostMap::performDownsampling(void* poDS_ptr, int full_w, int full_h, uint64_t bytes_needed, uint64_t target_bytes, const std::string& path) {
#ifdef HAVE_GDAL
    GDALDataset *poDS = (GDALDataset*)poDS_ptr;
    std::cerr << "ElevationMap: raster exceeds target size and no suitable overview found, performing in-code max downsample" << std::endl;
    GDALRasterBand *band = poDS->GetRasterBand(1);
    if (!band) {
      std::cerr << "ElevationMap: no raster band found" << std::endl;
      return false;
    }
    int hasNoData = 0;
    double noDataVal = band->GetNoDataValue(&hasNoData);

    double scale = std::sqrt((double)bytes_needed / (double)target_bytes);
    int factor = static_cast<int>(std::ceil(scale));
    if (factor < 1) factor = 1;

    const double MIN_VALID_FRACTION = 0.01;
    const int MAX_ITER = 8;
    int iter = 0;
    bool done = false;
    std::vector<float> best_data;
    int best_w = 0, best_h = 0;
    while (!done && iter < MAX_ITER) {
      ++iter;
      int out_w = static_cast<int>((full_w + factor - 1) / factor);
      int out_h = static_cast<int>((full_h + factor - 1) / factor);

      std::vector<float> outdata(static_cast<size_t>(out_w) * static_cast<size_t>(out_h), std::numeric_limits<float>::quiet_NaN());
      std::vector<char> hasValTmp(static_cast<size_t>(out_w) * static_cast<size_t>(out_h), 0);
      std::vector<float> rowbuf(full_w);
      for (int y = 0; y < full_h; ++y) {
        CPLErr rerr = band->RasterIO(GF_Read, 0, y, full_w, 1, rowbuf.data(), full_w, 1, GDT_Float32, 0, 0);
        if (rerr != CE_None) {
          std::cerr << "ElevationMap: RasterIO row read failed at y=" << y << "\n";
          return false;
        }
        int oy = y / factor;
        for (int x = 0; x < full_w; ++x) {
          float v = rowbuf[x];
          if (!std::isfinite(v)) continue;
          bool isNoData = false;
          if (hasNoData) {
            if (v == static_cast<float>(noDataVal)) isNoData = true;
          } else {
            if (v == -32767.0f || v == -32768.0f || v == -9999.0f || v == -99999.0f) isNoData = true;
          }
          if (isNoData) continue;
          int ox = x / factor;
          size_t idx = static_cast<size_t>(oy) * static_cast<size_t>(out_w) + static_cast<size_t>(ox);
          if (!hasValTmp[idx]) {
            outdata[idx] = v;
            hasValTmp[idx] = 1;
          } else {
            if (v > outdata[idx]) outdata[idx] = v;
          }
        }
      }
      size_t total_out = static_cast<size_t>(out_w) * static_cast<size_t>(out_h);
      size_t valid_count = 0;
      for (size_t i = 0; i < total_out; ++i) if (hasValTmp[i]) ++valid_count;
      double valid_frac = total_out ? (double)valid_count / (double)total_out : 0.0;

      if (valid_frac >= MIN_VALID_FRACTION || factor == 1) {
        best_data.swap(outdata);
        best_w = out_w; best_h = out_h;
        done = true;
        break;
      } else {
        best_data.swap(outdata);
        best_w = out_w; best_h = out_h;
        int new_factor = std::max(1, factor / 2);
        if (new_factor == factor) { done = true; break; }
        factor = new_factor;
      }
    }

    if (best_data.empty()) {
      std::cerr << "ElevationMap: downsample produced no data" << std::endl;
      return false;
    }

    elev_width_ = best_w; elev_height_ = best_h;
    elev_data_.swap(best_data);
    double eff_scale_x = static_cast<double>(full_w) / static_cast<double>(elev_width_);
    double eff_scale_y = static_cast<double>(full_h) / static_cast<double>(elev_height_);
    elev_geo_.pixel_w = elev_geo_.pixel_w * eff_scale_x;
    elev_geo_.pixel_h = elev_geo_.pixel_h * eff_scale_y;

    std::string ovrPath = path + ".ovr";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (poDriver) {
        char **papszOptions = nullptr;
        papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
        papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
        GDALDataset *poOvrDS = poDriver->Create(ovrPath.c_str(), elev_width_, elev_height_, 1, GDT_Float32, papszOptions);
        if (poOvrDS) {
             double new_gt[6];
             new_gt[0] = elev_geo_.origin_x;
             new_gt[1] = elev_geo_.pixel_w;
             new_gt[2] = elev_geo_.rot_x;
             new_gt[3] = elev_geo_.origin_y;
             new_gt[4] = elev_geo_.rot_y;
             new_gt[5] = elev_geo_.pixel_h;
             poOvrDS->SetGeoTransform(new_gt);
             poOvrDS->SetProjection(poDS->GetProjectionRef());
             GDALRasterBand *ovBand = poOvrDS->GetRasterBand(1);
             if (hasNoData) ovBand->SetNoDataValue(noDataVal);
             if (ovBand->RasterIO(GF_Write, 0, 0, elev_width_, elev_height_, elev_data_.data(), elev_width_, elev_height_, GDT_Float32, 0, 0) != CE_None) {
                 std::cerr << "ElevationMap: failed to write data to overview file " << ovrPath << std::endl;
             }
             GDALClose(poOvrDS);
        } else {
             std::cerr << "ElevationMap: failed to create overview file " << ovrPath << std::endl;
        }
        CSLDestroy(papszOptions);
    }

    elev_valid_ = true;
    std::cerr << "ElevationMap: in-code max downsample finished (width=" << elev_width_ << " height=" << elev_height_ << ")" << std::endl;
    return true;
#else
    return false;
#endif
}

void ElevationCostMap::printElevationInfo() const
{
  if (!elev_valid_) {
    std::cerr << "ElevationMap: no valid data loaded" << std::endl;
    return;
  }
  std::cerr << "ElevationMap: info -> width=" << elev_width_ << " height=" << elev_height_ << " resolution_x=" << elev_geo_.pixel_w << " resolution_y=" << elev_geo_.pixel_h << std::endl;
  double minv = std::numeric_limits<double>::infinity();
  double maxv = -std::numeric_limits<double>::infinity();
  uint64_t nan_count = 0;
  size_t total = static_cast<size_t>(elev_width_) * static_cast<size_t>(elev_height_);
  for (size_t i = 0; i < total; ++i) {
    float v = elev_data_[i];
    if (std::isnan(v)) { ++nan_count; continue; }
    if (v < minv) minv = v;
    if (v > maxv) maxv = v;
  }
  if (nan_count == total) {
    std::cerr << "ElevationMap: all values are NaN" << std::endl;
  } else {
    std::cerr << "ElevationMap: data min=" << minv << " max=" << maxv << " NaN_count=" << nan_count << std::endl;
  }
}

bool ElevationCostMap::getElevationAt(double x, double y, double &elev) const
{
  if (!elev_valid_) return false;
  double px = (x - elev_geo_.origin_x) / elev_geo_.pixel_w - 0.5;
  double py = (y - elev_geo_.origin_y) / elev_geo_.pixel_h - 0.5;
  int ix = static_cast<int>(floor(px));
  int iy = static_cast<int>(floor(py));
  double fx = px - ix;
  double fy = py - iy;
  if (ix < 0 || iy < 0 || ix+1 >= elev_width_ || iy+1 >= elev_height_) return false;
  int i00 = iy * elev_width_ + ix;
  int i10 = iy * elev_width_ + (ix+1);
  int i01 = (iy+1) * elev_width_ + ix;
  int i11 = (iy+1) * elev_width_ + (ix+1);
  double v00 = elev_data_[i00];
  double v10 = elev_data_[i10];
  double v01 = elev_data_[i01];
  double v11 = elev_data_[i11];
  elev = v00*(1-fx)*(1-fy) + v10*(fx)*(1-fy) + v01*(1-fx)*(fy) + v11*(fx)*(fy);
  return true;
}

void ElevationCostMap::createCostMap(int width, int height, double resolution, double origin_x, double origin_y)
{
  cost_width_ = width;
  cost_height_ = height;
  cost_resolution_ = resolution;
  cost_origin_x_ = origin_x;
  cost_origin_y_ = origin_y;
  cost_data_.assign(cost_width_ * cost_height_, -std::numeric_limits<float>::infinity());
}

void ElevationCostMap::setCost(int x, int y, float c)
{
  if (x < 0 || y < 0 || x >= cost_width_ || y >= cost_height_) return;
  cost_data_[y * cost_width_ + x] = c;
}

float ElevationCostMap::getCost(int x, int y) const
{
  if (x < 0 || y < 0 || x >= cost_width_ || y >= cost_height_) return -std::numeric_limits<float>::infinity();
  return cost_data_[y * cost_width_ + x];
}

bool ElevationCostMap::getCostAt(double x, double y, float &val) const
{
  int c = static_cast<int>(floor((x - cost_origin_x_) / cost_resolution_));
  int r = static_cast<int>(floor((cost_origin_y_ - y) / cost_resolution_));
  if (c < 0 || c >= cost_width_ || r < 0 || r >= cost_height_) return false;
  val = cost_data_[r * cost_width_ + c];
  return true;
}
