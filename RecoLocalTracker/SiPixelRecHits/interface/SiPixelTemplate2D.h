//
//  SiPixelTemplate2D.h (v1.03)
//
//  Full 2-D templates for cluster splitting, simulated cluster reweighting, and improved cluster probability
//
// Created by Morris Swartz on 12/01/09.
// V1.01 - fix qavg_ filling
// V1.02 - Add locBz to test if FPix use is out of range
// V1.03 - Fix edge checking on final template to increase template size and to properly truncate cluster
//
//

// Build the template storage structure from several pieces

#ifndef SiPixelTemplate2D_h
#define SiPixelTemplate2D_h 1

#include<vector>
#include<cassert>
#include "boost/multi_array.hpp"

#ifndef SI_PIXEL_TEMPLATE_STANDALONE
#include "CondFormats/SiPixelObjects/interface/SiPixelTemplateDBObject.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/SiPixelTemplateDefs.h"
#else
#include "SiPixelTemplateDefs.h"
#endif

struct SiPixelTemplateEntry2D { //!< Basic template entry corresponding to a single set of track angles
   int runnum;              //!< number of pixelav run used to generate this entry
   float cotalpha;          //!< cot(alpha) is proportional to cluster length in x and is basis of interpolation
   float cotbeta;           //!< cot(beta) is proportional to cluster length in y and is basis of interpolation
   float costrk[3];         //!< direction cosines of tracks used to generate this entry
   float qavg;              //!< average cluster charge for this set of track angles
   float pixmax;            //!< maximum charge for individual pixels in cluster
   float sxymax;            //!< average pixel signal for use of the error parameterization
   int iymin;               //!< the minimum nonzero pixel yindex in template (saves time during interpolation)
   int iymax;               //!< the maximum nonzero pixel yindex in template (saves time during interpolation)
   int jxmin;               //!< the minimum nonzero pixel xindex in template (saves time during interpolation)
   int jxmax;               //!< the maximum nonzero pixel xindex in template (saves time during interpolation)
   float xypar[2][5];       //!< pixel uncertainty parameterization
   float lanpar[2][5];      //!< pixel landau distribution parameters
   float xytemp[7][7][T2YSIZE][T2XSIZE];  //!< templates for y-reconstruction (binned over 1 central pixel)
   float chi2avg[4];        //!< average chi^2 in 4 charge bins
   float chi2min[4];        //!< minimum of chi^2 in 4 charge bins
   float chi2avgone;        //!< average y chi^2 for 1 pixel clusters
   float chi2minone;        //!< minimum of y chi^2 for 1 pixel clusters
   float clsleny;           //!< cluster y-length in pixels at signal height symax/2
   float clslenx;           //!< cluster x-length in pixels at signal height sxmax/2
   float spare[18];
} ;




struct SiPixelTemplateHeader2D {           //!< template header structure
   int ID;                 //!< template ID number
   int NTy;                //!< number of Template y entries
   int NTyx;               //!< number of Template y-slices of x entries
   int NTxx;               //!< number of Template x-entries in each slice
   int Dtype;              //!< detector type (0=BPix, 1=FPix)
   float qscale;           //!< Charge scaling to match cmssw and pixelav
   float lorywidth;        //!< estimate of y-lorentz width for optimal resolution
   float lorxwidth;        //!< estimate of x-lorentz width for optimal resolution
   float lorybias;         //!< estimate of y-lorentz bias
   float lorxbias;         //!< estimate of x-lorentz bias
   float Vbias;            //!< detector bias potential in Volts
   float temperature;      //!< detector temperature in deg K
   float fluence;          //!< radiation fluence in n_eq/cm^2
   float s50;              //!< 1/2 of the multihit dcol threshold in electrons
   float ss50;             //!< 1/2 of the single hit dcol threshold in electrons
   char title[80];         //!< template title
   int templ_version;      //!< Version number of the template to ensure code compatibility
   float Bfield;           //!< Bfield in Tesla
   float fbin[3];          //!< The QBin definitions in Q_clus/Q_avg
   float xsize;            //!< pixel size (for future use in upgraded geometry)
   float ysize;            //!< pixel size (for future use in upgraded geometry)
   float zsize;            //!< pixel size (for future use in upgraded geometry)
} ;



struct SiPixelTemplateStore2D { //!< template storage structure
   SiPixelTemplateHeader2D head;
#ifndef SI_PIXEL_TEMPLATE_USE_BOOST
   SiPixelTemplateEntry2D entry[61][5];  //!< use 2d entry to store [47][5] barrel entries or [5][9] fpix
#else
   boost::multi_array<SiPixelTemplateEntry2D,2> entry;  //!< use 2d entry to store [47][5] barrel entries or [5][9] fpix entries
#endif
} ;





// ******************************************************************************************
//! \class SiPixelTemplate2D
//!
//!  A template management class.  SiPixelTemplate contains thePixelTemp
//!  (a std::vector  of SiPixelTemplateStore, each of which is a collection of many
//!  SiPixelTemplateEntries).  Each SiPixelTemplateStore corresponds to a given detector
//!  condition, and is valid for a range of runs.  We allow more than one Store since the
//!  may change over time.
//!
//!  This class reads templates from files via pushfile() method.
//!
//!  The main functionality of SiPixelTemplate is xytemp(), which produces a template
//!  on the fly, given a specific track's alpha and beta.  The results are kept in data
//!  members and accessed via inline getters.
//!
//!  The resulting template is then used by PixelTempReco2D() (a global function) which
//!  get the reference for SiPixelTemplate & templ and uses the current template to
//!  reconstruct the SiPixelRecHit.
// ******************************************************************************************
class SiPixelTemplate2D {
public:
   SiPixelTemplate2D(const std::vector< SiPixelTemplateStore2D > & thePixelTemp) : thePixelTemp_(thePixelTemp) {id_current_ = -1; index_id_ = -1; cota_current_ = 0.; cotb_current_ = 0.;} //!< Default constructor
   
   static bool pushfile(int filenum, std::vector< SiPixelTemplateStore2D > & thePixelTemp_);     // load the private store with info from the
   // file with the index (int) filenum
   
#ifndef SI_PIXEL_TEMPLATE_STANDALONE
   static bool pushfile(const SiPixelTemplateDBObject& dbobject, std::vector< SiPixelTemplateStore2D > & thePixelTemp_);     // load the private store with info from db
#endif
   
   //  Initialize things before interpolating
   
   bool interpolate(int id, float cotalpha, float cotbeta, float locBz, float locBx);
   
   // Interpolate input alpha and beta angles to produce a working template for each individual hit.
   
   // Works with Phase 0+1
   bool xytemp(float xhit, float yhit, bool ydouble[BYM2], bool xdouble[BXM2], float template2d[BXM2][BYM2], bool dervatives, float dpdx2d[2][BXM2][BYM2]);
   
   // Overload for backward compatibility
   
   bool xytemp(float xhit, float yhit, bool ydouble[BYM2], bool xdouble[BXM2], float template2d[BXM2][BYM2]);
   
   void xysigma2(float qpixel, int index, float& xysig2);
   
   // Get the interpolated Landau distribution parameters
   
   void landau_par(float lanpar[2][5]);
   
   float qavg() {return qavg_;}        //!< average cluster charge for this set of track angles
   float pixmax() {return pixmax_;}    //!< maximum pixel charge
   float qscale() {return qscale_;}    //!< charge scaling factor
   float s50() {return s50_;}          //!< 1/2 of the pixel threshold signal in adc units
   float sxymax() {return sxymax_;}    //!< max pixel signal for pixel error calculation
   float chi2avg(int i) {
#ifndef SI_PIXEL_TEMPLATE_STANDALONE
      if(i < 0 || i > 3) {throw cms::Exception("DataCorrupt") << "SiPixelTemplate2D::chi2yavg called with illegal index = " << i << std::endl;}
#else
      assert(i>=0 && i<4);
#endif
      return chi2avg_[i];} //!< average chi^2 in 4 charge bins
   float chi2min(int i) {
#ifndef SI_PIXEL_TEMPLATE_STANDALONE
      if(i < 0 || i > 3) {throw cms::Exception("DataCorrupt") << "SiPixelTemplate2D::chi2ymin called with illegal index = " << i << std::endl;}
#else
      assert(i>=0 && i<4);
#endif
      return chi2min_[i];} //!< minimum chi^2 in 4 charge bins
   float fbin(int i) {
#ifndef SI_PIXEL_TEMPLATE_STANDALONE
      if(i < 0 || i > 2) {throw cms::Exception("DataCorrupt") << "SiPixelTemplate2D::fbin called with illegal index = " << i << std::endl;}
#else
      assert(i>=0 && i<3);
#endif
      return fbin_[i];} //!< Return lower bound of Qbin definition
   float sizex() {return clslenx_;}                             //! return x size of template cluster
   float sizey() {return clsleny_;}                             //! return y size of template cluster
   float chi2avgone() {return chi2avgone_;}                        //!< //!< average y chi^2 for 1 pixel clusters
   float chi2minone() {return chi2minone_;}                        //!< //!< minimum of y chi^2 for 1 pixel clusters
   float lorydrift() {return lorydrift_;}                            //!< signed lorentz y-width (microns)
   float lorxdrift() {return lorxdrift_;}                            //!< signed lorentz x-width (microns)
   float clsleny() {return clsleny_;}                                //!< cluster y-size
   float clslenx() {return clslenx_;}                                //!< cluster x-size
   float xsize() {return xsize_;}                                    //!< pixel x-size (microns)
   float ysize() {return ysize_;}                                    //!< pixel y-size (microns)
   float zsize() {return zsize_;}                                    //!< pixel z-size or thickness (microns)
   int storesize() {return (int)thePixelTemp_.size();}                    //!< return the size of the template store (the number of stored IDs
   
private:
   
   // Keep current template interpolaion parameters
   
   int id_current_;           //!< current id
   int index_id_;             //!< current index
   float cota_current_;       //!< current cot alpha
   float cotb_current_;       //!< current cot beta
   int Nyx_;                  //!< number of cot(beta)-entries (columns) in template
   int Nxx_;                  //!< number of cot(alpha)-entries (rows) in template
   int Dtype_;                //!< flags BPix (=0) or FPix (=1)
   float cotbeta0_;           //!< minimum cot(beta) covered
   float cotbeta1_;           //!< maximum cot(beta) covered
   float deltacotb_;          //!< cot(beta) bin size
   float cotalpha0_;          //!< minimum cot(alpha) covered
   float cotalpha1_;          //!< maximum cot(alpha) covered
   float deltacota_;          //!< cot(alpha) bin size
   int iy0_;                  //!< index of nearest cot(beta) bin
   int iy1_;                  //!< index of next-nearest cot(beta) bin
   float adcotb_;             //!< fractional pixel distance of cot(beta) from iy0_
   int jx0_;                  //!< index of nearest cot(alpha) bin
   int jx1_;                  //!< index of next-nearest cot(alpha) bin
   float adcota_;             //!< fractional pixel distance of cot(alpha) from jx0_
   int imin_;                 //!< min y index of templated cluster
   int imax_;                 //!< max y index of templated cluster
   int jmin_;                 //!< min x index of templated cluster
   int jmax_;                 //!< max x index of templated cluster
   bool flip_y_;              //!< flip y sign-sensitive quantities
   bool flip_x_;              //!< flip x sign-sensitive quantities
   bool success_;             //!< true if cotalpha, cotbeta are inside of the acceptance (dynamically loaded)
   
   
   // Keep results of last interpolation to return through member functions
   
   float qavg_;              //!< average cluster charge for this set of track angles
   float pixmax_;            //!< maximum pixel charge
   float qscale_;            //!< charge scaling factor
   float s50_;               //!< 1/2 of the pixel threshold signal in adc units
   float sxymax_;            //!< average pixel signal for y-projection of cluster
   float xytemp_[BXM2][BYM2];//!< template for xy-reconstruction
   float xypary0x0_[2][5];   //!< Polynomial error parameterization at ix0,iy0
   float xypary1x0_[2][5];   //!< Polynomial error parameterization at ix0,iy1
   float xypary0x1_[2][5];   //!< Polynomial error parameterization at ix1,iy0
   float lanpar_[2][5];      //!< Interpolated Landau parameters
   float chi2avg_[4];       //!< average chi^2 in 4 charge bins
   float chi2min_[4];       //!< minimum of chi^2 in 4 charge bins
   float chi2avgone_;       //!< average chi^2 for 1 pixel clusters
   float chi2minone_;       //!< minimum of chi^2 for 1 pixel clusters
   float clsleny_;            //!< projected y-length of cluster
   float clslenx_;            //!< projected x-length of cluster
   float lorywidth_;         //!< Lorentz y-width (sign corrected for fpix frame)
   float lorxwidth_;         //!< Lorentz x-width
   float lorydrift_;         //!< Lorentz y-drift
   float lorxdrift_;         //!< Lorentz x-drift
   float xsize_;             //!< Pixel x-size
   float ysize_;             //!< Pixel y-size
   float zsize_;             //!< Pixel z-size (thickness)
   float fbin_[3];           //!< The QBin definitions in Q_clus/Q_avg
   
   // The actual template store is a std::vector container
   
   const std::vector< SiPixelTemplateStore2D > & thePixelTemp_;
} ;


#endif
