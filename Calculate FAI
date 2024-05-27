var Lake_ROI = ee.FeatureCollection("projects/ee-haaaker4/assets/S_lakes");

//combine 250 and 500 resolution data
var addbands = function(image){

  var date = ee.Image(image).get('system:index');
  var time_start = ee.Image(image).get('time_start');
  var M250 = MOD250_rep.filter(ee.Filter.eq('system:index',date));
  
  var B1 = ee.Image(M250.get(0)).select(['sur_refl_b01'],['sur_refl_b01_250']);
  var B2 = ee.Image(M250.get(0)).select(['sur_refl_b02'],['sur_refl_b02_250']);
  var B3 = ee.Image(image).select('sur_refl_b03');


  image = ee.Image(image).addBands(B1);
  image = ee.Image(image).addBands(B2);

  image = ee.Image(image).updateMask(B3);
  
  return image.set('Date',date).set('time_start',time_start);
}

// helper function to extract the QA bits
function getQABits(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    // Return a single band image of the extracted QA bits, giving the band
    // a new name.
    return image.select([0], [newName])
                  .bitwiseAnd(pattern)
                  .rightShift(start);
}

// A function to mask out pixels that did not have observations.
var maskEmptyPixels = function(image) {
  // Find pixels that had observations.
  var withObs = image.select('num_observations_1km').gt(0)
  return image.updateMask(withObs)
}

// A function to mask out cloudy pixels.
var maskClouds = function(image) {
  // Select the QA band.
  var QA = image.select('state_1km')
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA,0, 2, 'internal_quality_flag');
  // Return an image masking out cloudy areas.
  return image.updateMask(internalQuality.eq(0));
}
// Add cloud mask using red band threshold
function cloudmask(image){
  var date = ee.Image(image).get('system:index');
  var time_start = ee.Image(image).get('time_start');
  var cloud = image.select('sur_refl_b01').lte(1500);
  return image.updateMask(cloud).set('Date',date).set('time_start',time_start);
}
//sensor zenith mask with flag
function zenithmask(image){
  var date = ee.Image(image).get('system:index');
  var time_start = ee.Image(image).get('time_start');
  var zenith = image.select('SensorZenith').lte(6000);
  return image.updateMask(zenith).set('Date',date).set('time_start',time_start);
}
//add area information
var fun_AddAreaPropery = function(image){
  var date = ee.Image(image).get('Date');
  var time_start = ee.Image(image).get('time_start');
  var Pixel_count = ee.Image(image).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: Feature,
    scale:250,
    maxPixels: 1e13
  });
  Pixel_count = ee.Number(Pixel_count.get('b06'));
  var area = Pixel_count.multiply(0.0625);
  var Lake_Area = ee.Number(Lake_area);
  var ratio = area.divide(Lake_Area).multiply(100);
  return ee.Image(image).set('Area_ratio',ratio).set('time_start',time_start)
};

//reprojection and resampling function of MOD
var fun_reproject = function(image){
    return ee.Image(image).setDefaultProjection('EPSG:4326',null,null);
  }

var fun_resample = function(image){
  return ee.Image(image).setDefaultProjection('EPSG:4326',null,250);
}
//calculate FAI image
var fun_FAI = function(image){
  var date = ee.Image(image).get('Date');
  var time_start = ee.Image(image).get('time_start');
  var BLUE_465 = image.select('BLUE_465');
  var GREEN_555 = image.select('GREEN_555');
  var RED_645 = image.select('RED_645');
  var NIR_859 = image.select('NIR_859');
  var SWIR_1240 = image.select('SWIR_1240');
  var SWIR_1640 = image.select('SWIR_1640');
  var SWIR_2130 = image.select('SWIR_2130');
  
 var NIR = image.expression('RED + (SWIR - RED) * ((859 - 645) / (1240 - 645))',{
  'RED':RED_645,
  'SWIR':SWIR_1240
});
  
  var FAI = NIR_859.subtract(NIR).rename('FAI');
  
  var MNDWI = image.expression('(GREEN-SWIR)/(GREEN+SWIR)',{
  'GREEN':GREEN_555,
  'SWIR':SWIR_1240
}).rename('MNDWI');
  var NDVI = image.expression('(RED-NIR)/(RED+NIR)',{
  'RED':RED_645,
  'NIR':NIR_859
}).rename('NDVI');
var CMI = image.expression('GREEN - BLUE - (SWIR - BLUE) * ((555 - 465) / (1240 - 465))',{
  'GREEN':GREEN_555,
  'BLUE':BLUE_465,
  'SWIR':SWIR_1240
}).rename('CMI');

  return FAI.set('Date',date).set('time_start',time_start).addBands([MNDWI,NDVI,CMI]);
}

//calculate unmask image
var fun_500_scale = function(image){
  var scale = 0.0001;
  var date = ee.Image(image).get('Date');
  var time_start = ee.Image(image).get('system:time_start');
  var BLUE_465 = ee.Image(image).select('sur_refl_b03').multiply(scale).rename('BLUE_465');
  var GREEN_555 = ee.Image(image).select('sur_refl_b04').multiply(scale).rename('GREEN_555');
  var RED_645 = ee.Image(image).select('sur_refl_b01').multiply(scale).rename('RED_645');
  var NIR_859 = ee.Image(image).select('sur_refl_b02').multiply(scale).rename('NIR_859');
  var SWIR_1240 = ee.Image(image).select('sur_refl_b05').multiply(scale).rename('SWIR_1240');
  var SWIR_1640 = ee.Image(image).select('sur_refl_b06').multiply(scale).rename('SWIR_1640');
  var SWIR_2130 = ee.Image(image).select('sur_refl_b07').multiply(scale).rename('SWIR_2130');
  
  // var BLUE_mask = BLUE_465.gte(0).mask(1);
  // var GREEN_mask = GREEN_555.gte(0).mask(1);
  var RED_mask = RED_645.gte(0).mask(1);
  var NIR_mask = NIR_859.gte(0).mask(1);
  var SWIR_mask = SWIR_1240.gte(0).mask(1);
  
  // var Blue_465 = BLUE_465.updateMask(BLUE_mask);
  // var Green_555 = GREEN_555.updateMask(GREEN_mask);
  var Red_645 = RED_645.updateMask(RED_mask);
  var Nir_859 = NIR_859.updateMask(NIR_mask);
  var Swir_1240 = SWIR_1240.updateMask(SWIR_mask);
  
  var image2 = ee.Image(BLUE_465).addBands(GREEN_555).addBands(Red_645).addBands(Nir_859).addBands(Swir_1240).addBands(SWIR_1640).addBands(SWIR_2130);
  
  
  return image2.set('Date',date).set('time_start',time_start);
}
//calculate unmask image
var fun_250_scale = function(image){
  var scale = 0.0001;
  var date = ee.Image(image).get('system:index');
  var time_start = ee.Image(image).get('system:time_start');
  var RED_645 = ee.Image(image).select('sur_refl_b01').multiply(scale).rename('RED_645_250');
  var NIR_859 = ee.Image(image).select('sur_refl_b02').multiply(scale).rename('NIR_859_250');
  
  var RED_mask = RED_645.gte(0).mask(1);
  var NIR_mask = NIR_859.gte(0).mask(1);

  var Red_645 = RED_645.updateMask(RED_mask);
  var Nir_859 = NIR_859.updateMask(NIR_mask);

  var image2 = ee.Image(Red_645).addBands(Nir_859);
  
  
  return image2.set('Date',date).set('time_start',time_start);
}

var fun_water = function(image){
  var MNDWI = image.select('MNDWI').rename('water');
  var NDVI = image.select('NDVI');
  var water = MNDWI.gte(NDVI).mask(1);
  return water;
}
//Apply LST scale
var fun_LST_scale = function(image){
  var scale = 0.02;
  var date = ee.Image(image).get('system:index');
  var time_start = ee.Image(image).get('system:time_start');
  var LST = ee.Image(image).select('LST_Day_1km').multiply(scale).rename('LST');
  
  LST = LST.subtract(273.15);
  
  return LST.set('Date',date).set('time_start',time_start);
}
//Add NDSI band
var addNDSI = function(image){
  var date = ee.Image(image).get('system:index');
  var time_start = ee.Image(image).get('time_start');
  var NDSI = ee.ImageCollection(NDSI_Index.filter(ee.Filter.eq('system:index', date))).first().select('NDSI');
  
  return image.set('Date',date).set('time_start',time_start).addBands(NDSI);
}
//mask Lake-Ice 
var icemask = function(image){
  var date = ee.Image(image).get('system:index');
  var time_start = ee.Image(image).get('time_start');
  var ice = ee.Image(image).select('NDSI').gte(0.4);
  // var B02 = ee.Image(ice).select('b02').gte(1100);
  // var B04 = ee.Image(B02).select('b04').gte(1000).updateMask(B02);
  var img = ee.Image(image).updateMask(ice);
  return img.set('Date',date).set('time_start',time_start);
}

/*------------------------------------Main function------------------------------------*/
/*------------------------------------Main function------------------------------------*/
/*------------------------------------Main function------------------------------------*/
/*------------------------------------Main function------------------------------------*/
/*------------------------------------Main function------------------------------------*/

var basinNames = Lake_ROI.reduceColumns(ee.Reducer.toList(), ["Lake_name"])
                       .get("list");
print(basinNames);

basinNames.evaluate(function(names){
for (var i=0; i < 1 ;i++){
var name = names[i];
print('name',name);
var feature = Lake_ROI.filter(ee.Filter.eq("Lake_name", name)).first();
var Feature = feature.geometry().buffer(-500);
var Lake_area = ee.Number(Feature.area().multiply(0.000001));
 
// Import MOD collection, change the date you want
var S_Year = '2000';
var S_Mon = '05';
var S_Day = '01';
var S_Date = ee.Date(S_Year+'-'+S_Mon+'-'+S_Day);

var E_Year = '2010';
var E_Mon = '12';
var E_Day = '31';
var E_Date = ee.Date(E_Year+'-'+E_Mon+'-'+E_Day);

   
    
var MOD250 = ee.ImageCollection('MODIS/061/MOD09GQ')
// .map(cloudmask)
// .map(maskEmptyPixels)
.filterBounds(Feature)
.filterDate(S_Date,E_Date);
    
var MOD500 = ee.ImageCollection('MODIS/061/MOD09GA')
.map(maskClouds)
.map(cloudmask)
.map(zenithmask)
.filterBounds(Feature)
.filterDate(S_Date,E_Date);
print('MOD500',MOD500);
var MOD500_unmask = ee.ImageCollection('MODIS/061/MOD09GA')
.map(maskClouds)
.map(cloudmask)
.map(zenithmask)
.filterBounds(Feature)
.filterDate(S_Date,E_Date);

var MOD10A1 = ee.ImageCollection('MODIS/061/MOD10A1')
.filterDate(S_Date,E_Date);
//.select('NDSI');
// print('MOD10A1',MOD10A1);
 
var MOD11A1 = ee.ImageCollection('MODIS/061/MOD11A1')
.filterDate(S_Date,E_Date)
.select('LST_Day_1km');
print('MOD11A1',MOD11A1);

var size_500 = ee.Number(MOD500.size());
var size_250 = ee.Number(MOD250.size());
var LST_Index = MOD11A1.map(fun_resample).map(fun_LST_scale).toList(size_500);
print('LST_Index',LST_Index);
//reproject and resample
var MOD250_rep = MOD250.toList(size_250);
var MOD500_res = MOD500.map(fun_resample).toList(size_500);
var MOD_500 = MOD500_res.map(fun_500_scale);
var MOD_250 = MOD250_rep.map(fun_250_scale);


var Terra_addLST = MOD_500.map(function(image){
  var date = ee.Image(image).get('Date');
  var time_start = ee.Image(image).get('time_start');
  
  var LST = ee.Image(ee.ImageCollection(LST_Index.filter(ee.Filter.eq('system:index',date))).first().select('LST'));
  var B1 = ee.Image(ee.ImageCollection(MOD_250.filter(ee.Filter.eq('Date', date))).first().select('RED_645_250'));
  var B2 = ee.Image(ee.ImageCollection(MOD_250.filter(ee.Filter.eq('Date', date))).first().select('NIR_859_250'));
  var img = ee.Image(image).addBands([LST,B1,B2]);
  return img.set('Date',date).set('time_start',time_start);
});
print('Terra_addLST',Terra_addLST);
var Terra_icemask = Terra_addLST.map(function(image){
  var date = ee.Image(image).get('Date');
  var time_start = ee.Image(image).get('time_start');
  var lst = ee.Image(image).select('LST').gte(0.4);

  var img = ee.Image(image);
  return img.set('Date',date).set('time_start',time_start);
});

var MOD_arearatio = Terra_icemask.map(function(image){
  var date = ee.Image(image).get('Date');
  var time_start = ee.Image(image).get('time_start');
  var Pixel_count = ee.Image(image).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: Feature,
    scale:250,
    maxPixels: 1e13
  });
  Pixel_count = ee.Number(Pixel_count.get('NIR_859'));
  var area = Pixel_count.multiply(0.0625);
  var Lake_Area = ee.Number(Lake_area);
  var ratio = area.divide(Lake_Area).multiply(100);
  return ee.Image(image).set('Area_ratio',ratio).set('time_start',time_start).set('Date',date);
});
print('MOD_arearatio',MOD_arearatio);
var Terra_Filter = MOD_arearatio.filter(ee.Filter.gte('Area_ratio',30));

    var rrs_num = Terra_Filter.size();

    //calculate FAI
    var MOD_FAI = Terra_Filter.map(function(image){
      
  var date = ee.Image(image).get('Date');
  var time_start = ee.Image(image).get('time_start');
  
  var BLUE_465 = ee.Image(image).select('BLUE_465');
  var GREEN_555 = ee.Image(image).select('GREEN_555');
  var RED_645 = ee.Image(image).select('RED_645_250');
  var NIR_859 = ee.Image(image).select('NIR_859_250');
  var SWIR_1240 = ee.Image(image).select('SWIR_1240');
  var SWIR_1640 = ee.Image(image).select('SWIR_1640');
  var SWIR_2130 = ee.Image(image).select('SWIR_2130');
  
  var NIR = ee.Image(image).expression('RED + (SWIR - RED) * ((859 - 645) / (1240 - 645))',{
  'RED':RED_645,
  'SWIR':SWIR_1240
});
  
  var FAI = NIR_859.subtract(NIR).rename('FAI');

  var img = ee.Image(image).addBands(FAI);
  
  return img.set('Date',date).set('time_start',time_start);
  });
  print('MOD_FAI',MOD_FAI);
  
  
 
}
})
