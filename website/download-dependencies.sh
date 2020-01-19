#!/usr/bin/env bash

cd css/lib

wget https://unpkg.com/tabulator-tables@4.2.3/dist/css/tabulator_simple.min.css
wget https://unpkg.com/tachyons@4.10.0/css/tachyons.min.css

cd ../..

cd scripts/lib

wget https://cdn.plot.ly/plotly-latest.min.js
wget https://cdn.jsdelivr.net/npm/vega@5.0.0/build/vega.js
wget https://cdn.jsdelivr.net/npm/vega-lite@3.0.0-rc14/build/vega-lite.js
wget https://cdn.jsdelivr.net/npm/vega-embed@3.29.1/build/vega-embed.js
wget https://unpkg.com/tabulator-tables@4.2.3/dist/js/tabulator.min.js
wget https://oss.sheetjs.com/js-xlsx/xlsx.full.min.js
wget https://d3js.org/d3.v5.min.js


