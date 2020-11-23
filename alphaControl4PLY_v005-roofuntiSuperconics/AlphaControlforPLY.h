#include <kvs/PointObject>
#ifndef _AlphaControlforPLY_H__
#define _AlphaControlforPLY_H__

#include <kvs/Module> 
#include <kvs/PointObject>
#include <kvs/PolygonObject> 

class AlphaControlforPLY: public kvs::PointObject {
  kvsModuleSuperClass( kvs::PointObject );

 public:
  AlphaControlforPLY( void );
  AlphaControlforPLY( kvs::PolygonObject* ply,
		      kvs::Camera* camera,
		      kvs::Vector3f BBMin,
		      kvs::Vector3f BBMax, 
		      int repeatLevel,
		      double alpha,
		      bool hface    );
  
 private:
  double m_searchRadius;
  double m_ratio;
  double m_pixel_width;
 private:
  void calculatePixelWidht( kvs::Camera* camera,
			    kvs::Vector3f BBMin, 
			    kvs::Vector3f BBMax );
  double calculateRequiredPartcleNumber( 
				      double  alpha, 
				      double repeatLevel,
				      kvs::Vector3f BBMin, 
				      kvs::Vector3f BBMax );
  void calculatePointRaio( const double analyticalNum, 
			   kvs::PolygonObject* ply );  
  void setParticles( kvs::PolygonObject* ply);


};

#endif
