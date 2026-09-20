// ==========================================
// MODERN NATURALIST
// Catalina Mountains Monsoon Vegetation
// ==========================================


// ------------------------------------------
// 1. Define our study area
// ------------------------------------------

var catalinas = ee.Geometry.Rectangle([
  -111.00, 32.20,
  -110.55, 32.60
]);


// ------------------------------------------
// 2. Get Sentinel-2 imagery
// ------------------------------------------

var image = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(catalinas)
  .filterDate('2026-08-01', '2026-09-01')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .median();


// ------------------------------------------
// 3. Calculate NDVI
// ------------------------------------------

var ndvi = image
  .normalizedDifference(['B8', 'B4'])
  .clip(catalinas);


// ------------------------------------------
// 4. Get elevation and aspect
// ------------------------------------------

var elevation = ee.Image('USGS/3DEP/10m');

var aspect = ee.Terrain.aspect(elevation);


// ------------------------------------------
// 5. Center map
// ------------------------------------------

Map.centerObject(catalinas, 11);


// ------------------------------------------
// 6. Display NDVI
// ------------------------------------------

Map.addLayer(
  ndvi,
  {
    min: 0,
    max: 0.8,
    palette: ['brown', 'yellow', 'green']
  },
  'NDVI'
);


// ------------------------------------------
// 7. Display elevation
// ------------------------------------------

Map.addLayer(
  elevation,
  {
    min: 500,
    max: 3000
  },
  'Elevation'
);


// ------------------------------------------
// 8. Display aspect
// ------------------------------------------

Map.addLayer(
  aspect,
  {
    min: 0,
    max: 360
  },
  'Aspect'
);


// ------------------------------------------
// 9. Define elevation bands
// ------------------------------------------

// Convert feet to meters

var lowMin = 3500 * 0.3048;
var lowMax = 5000 * 0.3048;

var midMin = 5000 * 0.3048;
var midMax = 7500 * 0.3048;

var highMin = 7500 * 0.3048;
var highMax = 9200 * 0.3048;


// Create elevation masks

var lowElevation = elevation
  .gte(lowMin)
  .and(elevation.lt(lowMax));

var midElevation = elevation
  .gte(midMin)
  .and(elevation.lt(midMax));

var highElevation = elevation
  .gte(highMin)
  .and(elevation.lte(highMax));


// ------------------------------------------
// 10. Define north and south slopes
// ------------------------------------------

// Aspect is meaningless on flat ground, so require a minimum slope
var slope = ee.Terrain.slope(elevation);
var steep = slope.gte(5);   // degrees; raise to 8 to be stricter

var north = aspect
  .gte(315)
  .or(aspect.lt(45))
  .and(steep);

var south = aspect
  .gte(135)
  .and(aspect.lt(225))
  .and(steep);

// ------------------------------------------
// 11. Create the six analysis buckets
// ------------------------------------------

var northLow = north.and(lowElevation);
var northMid = north.and(midElevation);
var northHigh = north.and(highElevation);

var southLow = south.and(lowElevation);
var southMid = south.and(midElevation);
var southHigh = south.and(highElevation);

// ------------------------------------------
// 12. Get June 9 Sentinel-2 imagery
// ------------------------------------------

var juneCollection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(catalinas)
  .filterDate('2026-06-09', '2026-06-10')
  .sort('CLOUDY_PIXEL_PERCENTAGE');

print('June 9 images:', juneCollection);
print('June 9 tiles:', juneCollection.aggregate_array('MGRS_TILE'));

// ------------------------------------------
// 13. Calculate June 9 NDVI
// ------------------------------------------

var juneImage = juneCollection.mosaic();

var juneNDVI = juneImage
  .normalizedDifference(['B8', 'B4'])
  .clip(catalinas);


// Display June 9 NDVI

Map.addLayer(
  juneNDVI,
  {
    min: 0,
    max: 0.8,
    palette: ['brown', 'yellow', 'green']
  },
  'June 9 NDVI'
);

// ------------------------------------------
// 14. Create zone image and zone mask
// ------------------------------------------

// 1 = North Low
// 2 = North Mid
// 3 = North High
// 4 = South Low
// 5 = South Mid
// 6 = South High

var zones = ee.Image(0)
  .where(northLow, 1)
  .where(northMid, 2)
  .where(northHigh, 3)
  .where(southLow, 4)
  .where(southMid, 5)
  .where(southHigh, 6)
  .rename('zone');

// Keep only pixels belonging to one of our six zones
var zoneMask = zones.gt(0);

// Mask the zone image
zones = zones.updateMask(zoneMask);

// ------------------------------------------
// 15. Cloud-masked NDVI, one image per acquisition
// ------------------------------------------

var START = '2026-06-01';
var WINDOW_DAYS = 14;
var N_WINDOWS = 7;     // 10 x 10 days = Jun 1 to Sep 9; keep within available data
var END = ee.Date(START).advance(WINDOW_DAYS * N_WINDOWS, 'day');

var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var QA_BAND = 'cs_cdf';
var CLEAR_THRESHOLD = 0.65;   // compositing makes a stricter mask affordable

var s2ndvi = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(catalinas)
  .filterDate(START, END)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
  .linkCollection(csPlus, [QA_BAND])
  .map(function(img) {
    var clear = img.select(QA_BAND).gte(CLEAR_THRESHOLD);
    return img
      .normalizedDifference(['B8', 'B4'])
      .rename('ndvi')
      .updateMask(clear)
      .copyProperties(img, ['system:time_start']);
  });

// ------------------------------------------
// 16. Median NDVI composite for each 10-day window
// ------------------------------------------

var startDate = ee.Date(START);

var windowStarts = ee.List.sequence(0, N_WINDOWS - 1).map(function(i) {
  return startDate.advance(ee.Number(i).multiply(WINDOW_DAYS), 'day');
});

// Label each window by its midpoint date (becomes the CSV column name)
var windowLabels = windowStarts.map(function(ws) {
  return ee.Date(ws).advance(WINDOW_DAYS / 2, 'day').format('YYYY-MM-dd');
});

var composites = ee.ImageCollection.fromImages(
  windowStarts.map(function(ws) {
    ws = ee.Date(ws);
    return s2ndvi
      .filterDate(ws, ws.advance(WINDOW_DAYS, 'day'))
      .median();
  })
);

var ndviStack = composites.toBands().rename(windowLabels);


// ------------------------------------------
// 17. Use the SAME pixels in every window
// ------------------------------------------

// Keep a pixel only if it has valid NDVI in all windows, so zone means
// change because vegetation changed, not because different pixels were visible.
var allValid = ndviStack.mask().reduce(ee.Reducer.min());
var analysisMask = allValid.and(zoneMask);
var ndviClean = ndviStack.updateMask(analysisMask);


// ------------------------------------------
// 21. Zone means per window, plus a pixel-retention check
// ------------------------------------------

var zoneResults = ee.FeatureCollection(
  ee.List.sequence(1, 6).map(function(z) {
    z = ee.Number(z);
    var stats = ndviClean
      .updateMask(zones.eq(z))
      .reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: catalinas,
        scale: 30,
        maxPixels: 1e13,
        tileScale: 8
      });
    return ee.Feature(null, stats).set('zone', z);
  })
);

// How many of each zone's pixels survived the "valid in every window" rule
var retention = ee.FeatureCollection(
  ee.List.sequence(1, 6).map(function(z) {
    z = ee.Number(z);
    var inZone = ee.Image(1).updateMask(zones.eq(z));
    var counts = inZone.rename('total')
      .addBands(inZone.updateMask(analysisMask).rename('kept'))
      .reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: catalinas,
        scale: 30,
        maxPixels: 1e13,
        tileScale: 8
      });
    return ee.Feature(null, counts).set('zone', z);
  })
);

Export.table.toDrive({
  collection: zoneResults,
  description: 'catalinas_ndvi_by_zone',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: retention,
  description: 'catalinas_ndvi_pixel_retention',
  fileFormat: 'CSV'
});
