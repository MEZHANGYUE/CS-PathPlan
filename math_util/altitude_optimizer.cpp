#include "altitude_optimizer.hpp"
#include <cmath>
#include <cpl_string.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <cpl_conv.h>
#endif

using namespace std;
namespace math_util {

ElevationMap::ElevationMap(): width_(0), height_(0), data_(), geo_(), valid_(false) {}
ElevationMap::~ElevationMap() {}

bool ElevationMap::loadFromTiff(const std::string &path)
{
#ifdef HAVE_GDAL
  GDALAllRegister();
  std::cerr << "ElevationMap: opening TIFF: " << path << std::endl;
  // open readonly first to probe size and geotransform
  GDALDataset *poDS = (GDALDataset*)GDALOpen(path.c_str(), GA_ReadOnly);
  if (!poDS) {
    cerr << "ElevationMap: GDALOpen failed for " << path << "\n";
    return false;
  }
  int full_w = poDS->GetRasterXSize();
  int full_h = poDS->GetRasterYSize();
  std::cerr << "ElevationMap: dataset size: " << full_w << " x " << full_h << std::endl;
  double gt[6];
  if (poDS->GetGeoTransform(gt) == CE_None) {
    geo_.origin_x = gt[0];
    geo_.pixel_w = gt[1];
    geo_.rot_x = gt[2];
    geo_.origin_y = gt[3];
    geo_.rot_y = gt[4];
    geo_.pixel_h = gt[5];
    std::cerr << "ElevationMap: geotransform: origin(" << geo_.origin_x << ", " << geo_.origin_y << ") pixel_w=" << geo_.pixel_w << " pixel_h=" << geo_.pixel_h << std::endl;
  } else {
    // default geotransform (pixel-as-meter)
    geo_.origin_x = 0.0; geo_.pixel_w = 1.0; geo_.rot_x = 0.0;
    geo_.origin_y = 0.0; geo_.rot_y = 0.0; geo_.pixel_h = -1.0;
  }

  // estimate memory needed (float32 per pixel)
  const uint64_t bytes_needed = static_cast<uint64_t>(full_w) * static_cast<uint64_t>(full_h) * sizeof(float);
   // target maximum bytes after downsampling: 1 GB (user request)
   const uint64_t TARGET_BYTES = 1ULL * 1024ULL * 1024ULL * 1024ULL;
   std::cerr << "ElevationMap: estimated memory for full raster (bytes): " << bytes_needed << " target_bytes=" << TARGET_BYTES << std::endl;

  // Check for existing overviews first
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

   // If file is large and no suitable overview exists, perform in-code downsampling using MAX pooling
   if (bytes_needed > TARGET_BYTES && !found_suitable_ov) {
     bool ok = performDownsampling(poDS, full_w, full_h, bytes_needed, TARGET_BYTES, path);
     GDALClose(poDS);
     if (ok) {
       printInfo();
       return true;
     } else {
       return false;
     }
   }  // If file is large and no overviews, try to create external overviews (.ovr) by reopening in update mode
  if (bytes_needed > TARGET_BYTES && existing_ov == 0) {
    std::cerr << "ElevationMap: large raster and no overviews; attempting to build external overviews..." << std::endl;
  GDALClose(poDS);
  // try to open in update mode to build overviews
  // Force creation of external overviews (.ovr) by setting GDAL config option
  CPLSetConfigOption("GDAL_TIFF_OVR", "YES");
  std::cerr << "ElevationMap: CPLSetConfigOption GDAL_TIFF_OVR=YES to force external overviews" << std::endl;
  GDALDataset *poDSup = (GDALDataset*)GDALOpen(path.c_str(), GA_Update);
    if (poDSup) {
      // build overview levels 2,4,8,... until smallest dimension < 256
      std::vector<int> ov_levels;
      int w = full_w, h = full_h;
      int level = 2;
      while ((w / level) >= 25600 || (h / level) >= 25600) {
        ov_levels.push_back(level);
        level *= 2;
      }
      std::cerr << "ElevationMap: requested overview levels: ";
      for (size_t i=0;i<ov_levels.size();++i) std::cerr << ov_levels[i] << (i+1<ov_levels.size()?",":"\n");
  if (!ov_levels.empty()) {
        // BuildOverviews signature expects: resampling, nListCount, panOverviewList, nBands, panBands, pfnProgress, pProgressData
        int nBands = 1;
        int panBands[1] = {1};
        // Try to build overviews using MAX sampling (take maximum of covered pixels) to preserve peaks.
        const char *resampling_try = "MAX";
        CPLErr berr = poDSup->BuildOverviews(resampling_try, static_cast<int>(ov_levels.size()), ov_levels.data(), nBands, panBands, nullptr, nullptr);
        if (berr != CE_None) {
          std::cerr << "ElevationMap: BuildOverviews with resampling=MAX failed, falling back to AVERAGE for " << path << "\n";
          resampling_try = "AVERAGE";
          berr = poDSup->BuildOverviews(resampling_try, static_cast<int>(ov_levels.size()), ov_levels.data(), nBands, panBands, nullptr, nullptr);
        }
        if (berr != CE_None) {
          cerr << "ElevationMap: BuildOverviews failed for " << path << "\n";
        } else {
          std::cerr << "ElevationMap: external overviews built for " << path << " using resampling=" << resampling_try << "\n";
        }
      } else {
        std::cerr << "ElevationMap: no overview levels chosen (small raster)" << std::endl;
      }
      GDALClose(poDSup);
    } else {
      // cannot open update; warn and fall back to read-only
      std::cerr << "ElevationMap: cannot open " << path << " in update mode to build overviews; continuing" << std::endl;
    }
    // clear the config option we set earlier
    CPLSetConfigOption("GDAL_TIFF_OVR", nullptr);
    // reopen readonly to continue
    poDS = (GDALDataset*)GDALOpen(path.c_str(), GA_ReadOnly);
    if (!poDS) {
      cerr << "ElevationMap: GDALOpen failed for " << path << " after overview attempt\n";
      return false;
    }
    band = poDS->GetRasterBand(1);
    existing_ov = band ? band->GetOverviewCount() : 0;
  }

  // If overviews exist, try to pick a lower-resolution overview that fits memory threshold
  GDALRasterBand *readBand = band;
  int use_ov_index = -1;
  if (band && existing_ov > 0) {
    // iterate from largest overview (highest index) to smallest
    for (int i = existing_ov - 1; i >= 0; --i) {
      GDALRasterBand *ob = band->GetOverview(i);
      if (!ob) continue;
      uint64_t ow = static_cast<uint64_t>(ob->GetXSize());
      uint64_t oh = static_cast<uint64_t>(ob->GetYSize());
  uint64_t need = ow * oh * sizeof(float);
  if (need <= TARGET_BYTES) {
        readBand = ob;
        use_ov_index = i;
        width_ = static_cast<int>(ow);
        height_ = static_cast<int>(oh);
        // adjust geotransform: pixel size scales by factor = full_w / ow (assume integer)
        double scale_x = static_cast<double>(full_w) / static_cast<double>(ow);
        double scale_y = static_cast<double>(full_h) / static_cast<double>(oh);
        geo_.pixel_w = geo_.pixel_w * scale_x;
        geo_.pixel_h = geo_.pixel_h * scale_y;
        break;
      }
    }
    std::cerr << "ElevationMap: existing overviews count=" << existing_ov << " chosen overview index=" << use_ov_index << std::endl;
  }

  // If we didn't pick an overview, read full resolution
  if (use_ov_index == -1) {
    width_ = full_w;
    height_ = full_h;
  }

  // allocate and read from the chosen band (overview or base)
  data_.assign(static_cast<size_t>(width_) * static_cast<size_t>(height_), std::numeric_limits<float>::quiet_NaN());
  CPLErr err;
  if (readBand) {
    // readBand may be an overview; use RasterIO on the band
    err = readBand->RasterIO(GF_Read, 0, 0, width_, height_, data_.data(), width_, height_, GDT_Float32, 0, 0);
  } else {
    err = CE_Failure;
  }
  if (err != CE_None) {
    cerr << "ElevationMap: GDAL RasterIO failed\n";
    GDALClose(poDS);
    return false;
  }
  GDALClose(poDS);
  valid_ = true;
  std::cerr << "ElevationMap: finished reading raster (width=" << width_ << " height=" << height_ << ")" << std::endl;
  // print a short summary of the map
  printInfo();
  return true;
#else
  cerr << "ElevationMap: GDAL not available at build time. Cannot load TIFF.\n";
  return false;
#endif
}

bool ElevationMap::performDownsampling(void* poDS_ptr, int full_w, int full_h, uint64_t bytes_needed, uint64_t target_bytes, const std::string& path) {
#ifdef HAVE_GDAL
    GDALDataset *poDS = (GDALDataset*)poDS_ptr;
    std::cerr << "ElevationMap: raster exceeds target size and no suitable overview found, performing in-code max downsample\n";
    GDALRasterBand *band = poDS->GetRasterBand(1);
    if (!band) {
      std::cerr << "ElevationMap: no raster band found\n";
      return false;
    }
    int hasNoData = 0;
    double noDataVal = band->GetNoDataValue(&hasNoData);

    // compute initial reduction factor so that output fits within TARGET_BYTES
    double scale = std::sqrt((double)bytes_needed / (double)target_bytes);
    int factor = static_cast<int>(std::ceil(scale));
    if (factor < 1) factor = 1;

    // We'll try decreasing the factor if the resulting downsample has too little valid data
    const double MIN_VALID_FRACTION = 0.01; // require at least 1% valid pixels
    const int MAX_ITER = 8;
    int iter = 0;
    bool done = false;
    std::vector<float> best_data;
    int best_w = 0, best_h = 0;
    while (!done && iter < MAX_ITER) {
      ++iter;
      int out_w = static_cast<int>((full_w + factor - 1) / factor);
      int out_h = static_cast<int>((full_h + factor - 1) / factor);
      std::cerr << "ElevationMap: attempt " << iter << " downsample factor=" << factor << " out_size=" << out_w << "x" << out_h << std::endl;

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
      std::cerr << "ElevationMap: valid fraction=" << valid_frac << " (" << valid_count << "/" << total_out << ")" << std::endl;

      if (valid_frac >= MIN_VALID_FRACTION || factor == 1) {
        // accept this result
        best_data.swap(outdata);
        best_w = out_w; best_h = out_h;
        done = true;
        break;
      } else {
        // keep as fallback, but try to reduce factor to get more valid pixels
        best_data.swap(outdata);
        best_w = out_w; best_h = out_h;
        int new_factor = std::max(1, factor / 2);
        if (new_factor == factor) { done = true; break; }
        factor = new_factor;
        std::cerr << "ElevationMap: too few valid pixels, reducing factor to " << factor << " and retrying" << std::endl;
      }
    }

    if (best_data.empty()) {
      std::cerr << "ElevationMap: downsample produced no data" << std::endl;
      return false;
    }

    // move best_data into data_
    width_ = best_w; height_ = best_h;
    data_.swap(best_data);
    // update geo transform based on final factor estimate: compute effective factor = full_w / width_
    double eff_scale_x = static_cast<double>(full_w) / static_cast<double>(width_);
    double eff_scale_y = static_cast<double>(full_h) / static_cast<double>(height_);
    geo_.pixel_w = geo_.pixel_w * eff_scale_x;
    geo_.pixel_h = geo_.pixel_h * eff_scale_y;

    // Save the downsampled data to .ovr file for future use
    std::string ovrPath = path + ".ovr";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (poDriver) {
        char **papszOptions = nullptr;
        papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
        papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
        GDALDataset *poOvrDS = poDriver->Create(ovrPath.c_str(), width_, height_, 1, GDT_Float32, papszOptions);
        if (poOvrDS) {
             double new_gt[6];
             new_gt[0] = geo_.origin_x;
             new_gt[1] = geo_.pixel_w;
             new_gt[2] = geo_.rot_x;
             new_gt[3] = geo_.origin_y;
             new_gt[4] = geo_.rot_y;
             new_gt[5] = geo_.pixel_h;
             poOvrDS->SetGeoTransform(new_gt);
             poOvrDS->SetProjection(poDS->GetProjectionRef());
             GDALRasterBand *ovBand = poOvrDS->GetRasterBand(1);
             if (hasNoData) ovBand->SetNoDataValue(noDataVal);
             if (ovBand->RasterIO(GF_Write, 0, 0, width_, height_, data_.data(), width_, height_, GDT_Float32, 0, 0) != CE_None) {
                 std::cerr << "ElevationMap: failed to write data to overview file " << ovrPath << std::endl;
             }
             GDALClose(poOvrDS);
             std::cerr << "ElevationMap: saved downsampled data to " << ovrPath << std::endl;
        } else {
             std::cerr << "ElevationMap: failed to create overview file " << ovrPath << std::endl;
        }
        CSLDestroy(papszOptions);
    }

    valid_ = true;
    std::cerr << "ElevationMap: in-code max downsample finished (width=" << width_ << " height=" << height_ << ")" << std::endl;
    return true;
#else
    return false;
#endif
}

void ElevationMap::printInfo() const
{
  if (!valid_) {
    std::cerr << "ElevationMap: no valid data loaded" << std::endl;
    return;
  }
  std::cerr << "ElevationMap: info -> width=" << width_ << " height=" << height_ << " resolution_x=" << geo_.pixel_w << " resolution_y=" << geo_.pixel_h << std::endl;
  // compute basic stats (skip NaNs)
  double minv = std::numeric_limits<double>::infinity();
  double maxv = -std::numeric_limits<double>::infinity();
  uint64_t nan_count = 0;
  size_t total = static_cast<size_t>(width_) * static_cast<size_t>(height_);
  for (size_t i = 0; i < total; ++i) {
    float v = data_[i];
    if (std::isnan(v)) { ++nan_count; continue; }
    if (v < minv) minv = v;
    if (v > maxv) maxv = v;
  }
  if (nan_count == total) {
    std::cerr << "ElevationMap: all values are NaN" << std::endl;
  } else {
    std::cerr << "ElevationMap: data min=" << minv << " max=" << maxv << " NaN_count=" << nan_count << std::endl;
  }
  // print a few sample values
  std::cerr << "ElevationMap: sample values (first up to 8): ";
  for (size_t i = 0; i < std::min<size_t>(8, total); ++i) {
    if (i) std::cerr << ", ";
    std::cerr << data_[i];
  }
  std::cerr << std::endl;
}

bool ElevationMap::loadFromPGM(const std::string &path)
{
  ifstream ifs(path, ios::binary);
  if (!ifs.is_open()) return false;
  string magic;
  ifs >> magic;
  if (magic != "P5" && magic != "P2") {
    cerr << "ElevationMap: Unsupported PGM format: " << magic << "\n";
    return false;
  }
  int w,h,maxv;
  ifs >> w >> h >> maxv;
  width_ = w; height_ = h;
  data_.assign(width_ * height_, 0.0f);
  if (magic == "P5") {
    ifs.get(); // consume single whitespace
    vector<unsigned char> buf(width_ * height_);
    ifs.read((char*)buf.data(), buf.size());
    for (int i=0;i<width_*height_;++i) data_[i] = static_cast<float>(buf[i]);
  } else {
    for (int i=0;i<width_*height_;++i) { int v; ifs >> v; data_[i] = static_cast<float>(v); }
  }
  geo_.origin_x = 0; geo_.origin_y = 0; geo_.pixel_w = 1.0; geo_.pixel_h = -1.0;
  valid_ = true;
  return true;
}

bool ElevationMap::getElevationAt(double x, double y, double &elev) const
{
  if (!valid_) return false;
  // Convert world x,y to pixel indices (assuming geo_.origin is top-left)
  // pixel center at origin + (i+0.5)*pixel_w
  double px = (x - geo_.origin_x) / geo_.pixel_w - 0.5;
  double py = (y - geo_.origin_y) / geo_.pixel_h - 0.5; // pixel_h may be negative
  // Note: if pixel_h negative, dividing by negative flips sign; above handles both
  // Bilinear interpolation
  int ix = static_cast<int>(floor(px));
  int iy = static_cast<int>(floor(py));
  double fx = px - ix;
  double fy = py - iy;
  if (ix < 0 || iy < 0 || ix+1 >= width_ || iy+1 >= height_) return false;
  int i00 = iy * width_ + ix;
  int i10 = iy * width_ + (ix+1);
  int i01 = (iy+1) * width_ + ix;
  int i11 = (iy+1) * width_ + (ix+1);
  double v00 = data_[i00];
  double v10 = data_[i10];
  double v01 = data_[i01];
  double v11 = data_[i11];
  elev = v00*(1-fx)*(1-fy) + v10*(fx)*(1-fy) + v01*(1-fx)*(fy) + v11*(fx)*(fy);
  return true;
}

// CostMap implementation
void CostMap::create(int width, int height, double resolution, double origin_x, double origin_y)
{
  width_ = width; height_ = height; resolution_ = resolution; origin_x_ = origin_x; origin_y_ = origin_y;
  costs_.assign(width_ * height_, 255);
}

void CostMap::setCost(int x, int y, unsigned char c)
{
  if (x < 0 || y < 0 || x >= width_ || y >= height_) return;
  costs_[y * width_ + x] = c;
}

unsigned char CostMap::getCost(int x, int y) const
{
  if (x < 0 || y < 0 || x >= width_ || y >= height_) return 255;
  return costs_[y * width_ + x];
}

void CostMap::fromElevationMap(const ElevationMap &emap, double elev_min, double elev_max)
{
  if (!emap.isValid()) return;
  int w = emap.getWidth();
  int h = emap.getHeight();
  double res = emap.getResolution();
  GeoTransform gt = emap.getGeoTransform();
  create(w, h, res, gt.origin_x, gt.origin_y + gt.pixel_h * h); // origin set to bottom-left
  for (int r = 0; r < h; ++r) {
    for (int c = 0; c < w; ++c) {
      int idx = r * w + c;
      double elev;
      // compute world coords for pixel center
      double wx = gt.origin_x + (c + 0.5) * gt.pixel_w;
      double wy = gt.origin_y + (r + 0.5) * gt.pixel_h;
      if (!emap.getElevationAt(wx, wy, elev)) {
        setCost(c, r, 255);
        continue;
      }
      double scaled = 0.0;
      if (elev_max > elev_min) scaled = (elev - elev_min) / (elev_max - elev_min);
      int val = static_cast<int>(round(scaled * 254.0));
      if (val < 0) val = 0; if (val > 254) val = 254;
      setCost(c, r, static_cast<unsigned char>(val));
    }
  }
}

// AltitudeOptimizer implementation
AltitudeOptimizer::AltitudeOptimizer() {}

bool AltitudeOptimizer::load_elevation_data(const std::string &elevation_file)
{
  // choose loader by extension
  std::string lower = elevation_file;
  for (auto &c : lower) c = static_cast<char>(std::tolower(c));
  bool ok = false;
  if (lower.size() >= 4 && (lower.rfind(".tif") == lower.size() - 4 || lower.rfind(".tiff") == lower.size() - 5)) {
    ok = elevation_map_.loadFromTiff(elevation_file);
    if (ok) {
      std::cout << "Elevation data loaded from: " << elevation_file << std::endl;
      return true;
    }
  }
  // try PGM
  if (lower.size() >= 4 && (lower.rfind(".pgm") == lower.size() - 4)) {
    ok = elevation_map_.loadFromPGM(elevation_file);
    if (ok) {
      std::cout << "Elevation data loaded from: " << elevation_file << std::endl;
      return true;
    }
  }
  // fallback: try both loaders
  if (!ok) {
    ok = elevation_map_.loadFromTiff(elevation_file);
    if (ok) {
      std::cout << "Elevation data loaded from: " << elevation_file << std::endl;
      return true;
    }
    ok = elevation_map_.loadFromPGM(elevation_file);
    if (ok) {
      std::cout << "Elevation data loaded from: " << elevation_file << std::endl;
      return true;
    }
  }
  if (!ok) {
    std::cerr << "Elevation data load failed for: " << elevation_file << std::endl;
  }
  return ok;
}

bool AltitudeOptimizer::optimizeHeights(const std::vector<Eigen::Vector3d> &waypoints,
                                       const ElevationMap &emap,
                                       const CostMap &cmap,
                                       const Params &p,
                                       std::vector<double> &out_z) const
{
  size_t n = waypoints.size();
  if (n == 0) return false;
  // Build Hessian H and rhs b for quadratic objective: 0.5*z^T H z - b^T z
  Eigen::SparseMatrix<double> H(n, n);
  std::vector<Eigen::Triplet<double>> trips;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(n);

  // smoothing: discrete second derivative penalty -> use operator L where (L z)_i = z_{i-1} - 2 z_i + z_{i+1}
  // For interior points i=1..n-2, add lambda_smooth * (L^T L) contributions
  if (n >= 3 && p.lambda_smooth > 0.0) {
    for (size_t i = 1; i + 1 < n; ++i) {
      // entries for L^T L at rows/cols i-1,i,i+1
      // contribution: [1 -2 1]^T * [1 -2 1] scaled by lambda
      // which yields a 3x3 block:
      // [1 -2 1]
      // [-2 4 -2]
      // [1 -2 1]
      double s = p.lambda_smooth;
      trips.emplace_back(i-1, i-1, s * 1.0);
      trips.emplace_back(i-1, i,   s * -2.0);
      trips.emplace_back(i-1, i+1, s * 1.0);

      trips.emplace_back(i, i-1,   s * -2.0);
      trips.emplace_back(i, i,     s * 4.0);
      trips.emplace_back(i, i+1,   s * -2.0);

      trips.emplace_back(i+1, i-1, s * 1.0);
      trips.emplace_back(i+1, i,   s * -2.0);
      trips.emplace_back(i+1, i+1, s * 1.0);
    }
  }

  // terrain follow weight: lambda_follow * (z - elev)^2 -> contributes lambda_follow * I to H and lambda_follow * elev to b
  for (size_t i = 0; i < n; ++i) {
    double wx = waypoints[i](0);
    double wy = waypoints[i](1);
    double elev = 0.0;
    if (emap.getElevationAt(wx, wy, elev)) {
      double s = p.lambda_follow;
      // follow target is terrain + UAV clearance (uav_R)
      double target = elev + p.uav_R;
      trips.emplace_back(i, i, s);
      b(i) += s * target;
    }
  }

  // max climb rate soft-constraint: penalize large vertical differences between consecutive points
  // term for edge (i,i+1): ((z_{i+1}-z_i)/(dist * max_climb_rate))^2
  // which expands to weight * (z_i^2 - 2 z_i z_{i+1} + z_{i+1}^2)
  if (p.max_climb_rate > 0.0) {
    for (size_t i = 0; i + 1 < n; ++i) {
      double x0 = waypoints[i](0), y0 = waypoints[i](1);
      double x1 = waypoints[i+1](0), y1 = waypoints[i+1](1);
      double dist = std::hypot(x1 - x0, y1 - y0);
      if (dist <= 1e-9) continue;
      double denom = (dist * p.max_climb_rate);
      if (denom <= 1e-12) continue;
      double w = 1.0 / (denom * denom); // weight for (z_{i+1}-z_i)^2
      trips.emplace_back(i, i,     w);
      trips.emplace_back(i, i+1,  -w);
      trips.emplace_back(i+1, i,  -w);
      trips.emplace_back(i+1, i+1, w);
    }
  }

  // add tiny diagonal regularization triplets
  for (size_t i = 0; i < n; ++i) {
    trips.emplace_back(i, i, 1e-8);
  }
  // assemble H
  H.setFromTriplets(trips.begin(), trips.end());

  // Solve H z = b
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(H);
  if (solver.info() != Eigen::Success) {
    cerr << "AltitudeOptimizer: solver decomposition failed" << endl;
    return false;
  }
  Eigen::VectorXd z = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    cerr << "AltitudeOptimizer: solver failed" << endl;
    return false;
  }
  out_z.resize(n);
  for (size_t i = 0; i < n; ++i) out_z[i] = z(i);
  return true;
}

} // namespace math_util
