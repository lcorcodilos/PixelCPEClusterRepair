#ifndef RecoLocalTracker_SiPixelRecHits_PixelCPEClusterRepair_H
#define RecoLocalTracker_SiPixelRecHits_PixelCPEClusterRepair_H

#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEBase.h"

// Already in the base class
//#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
//#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
//#include "Geometry/CommonDetAlgo/interface/MeasurementPoint.h"
//#include "Geometry/CommonDetAlgo/interface/MeasurementError.h"
//#include "Geometry/Surface/interface/GloballyPositioned.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"


#ifndef SI_PIXEL_TEMPLATE_STANDALONE
#include "RecoLocalTracker/SiPixelRecHits/interface/SiPixelTemplate2D.h"
#else
#include "SiPixelTemplate2D.h"
#endif

#include <utility>
#include <vector>


#if 0
/** \class PixelCPEClusterRepair
 * Perform the position and error evaluation of pixel hits using
 * 2D templates and local angles, via a 2D template morphing fit.
 */
#endif

class MagneticField;
class PixelCPEClusterRepair : public PixelCPEBase
{
public:
   struct ClusterParamTemplate2D : ClusterParam
   {
      ClusterParamTemplate2D(const SiPixelCluster & cl) : ClusterParam(cl){}
      // The result of PixelTemplateReco2D
      float templXrec_ ;
      float templYrec_ ;
      float templSigmaX_ ;
      float templSigmaY_ ;
      // Add new information produced by SiPixelTemplateReco::PixelTempReco2D &&&
      // These can only be accessed if we change silicon pixel data formats and add them to the rechit
      float templProbXY_ ;
      
      float templProbQ_;
      
      int templQbin_ ;
      
      int ierr;
      
   };
   
   // PixelCPEClusterRepair( const DetUnit& det );
   PixelCPEClusterRepair(edm::ParameterSet const& conf, const MagneticField *, const TrackerGeometry&, const TrackerTopology&,
                        const SiPixelLorentzAngle *, const SiPixelTemplateDBObject *);
   
   ~PixelCPEClusterRepair();
   
private:
   ClusterParam * createClusterParam(const SiPixelCluster & cl) const;
   
   // Compute the local position
   LocalPoint localPosition (DetParam const & theDetParam, ClusterParam & theClusterParam) const;
   
   // Compute the local errors
   LocalError localError   (DetParam const & theDetParam, ClusterParam & theClusterParam) const;
   
   // Template storage
   std::vector< SiPixelTemplateStore2D > thePixelTemp_;
   
   int speed_ ;
   
   bool UseClusterSplitter_;
   
   //bool DoCosmics_;
   //bool LoadTemplatesFromDB_;
   
};

#endif




