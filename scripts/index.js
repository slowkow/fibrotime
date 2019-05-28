
var fetch = function(data_url, callback) {
  logTime('downloading ' + data_url);
  var oReq = new XMLHttpRequest();
  oReq.open("GET", data_url, true);
  oReq.responseType = "arraybuffer";
  oReq.onload = function(oEvent) {
    var arrayBuffer = oReq.response;
    if (arrayBuffer) {
      var byteArray = new Uint8Array(arrayBuffer); // Get the compressed data
      byteArray = pako.inflate(byteArray) // Decompress the data
      var s = new TextDecoder("utf-8").decode(byteArray) // Convert to string
      callback(s)
    }
  };
  oReq.send(null);
}

var logTime = function(text) {
  console.log('[' + new Date().toUTCString() + '] ' + text);
}

// VEGA FIGURE 

// Assign the specification to a local variable vlSpec.
//var vlSpec = {
//  $schema: 'https://vega.github.io/schema/vega-lite/v3.json',
//  data: {
//    values: [
//      {a: 'C', b: 2},
//      {a: 'C', b: 7},
//      {a: 'C', b: 4},
//      {a: 'D', b: 1},
//      {a: 'D', b: 2},
//      {a: 'D', b: 6},
//      {a: 'E', b: 8},
//      {a: 'E', b: 4},
//      {a: 'E', b: 7}
//    ]
//  },
//  mark: 'bar',
//  encoding: {
//    y: {field: 'a', type: 'nominal'},
//    x: {
//      aggregate: 'average',
//      field: 'b',
//      type: 'quantitative',
//      axis: {
//        title: 'Average of b'
//      }
//    }
//  }
//};

// Assign the specification to a local variable vlSpec.
var vlSpec = {
  // $schema: 'https://vega.github.io/schema/vega-lite/v3.json',
  $schema: 'https://vega.github.io/schema/vega-lite/v3.0.0-rc14.json',
  data: {
    values: mydata,
  },
  mark: 'bar',
  encoding: {
    y: {field: 'stimulation', type: 'nominal'},
    x: {
      aggregate: 'average',
      field: 'gene',
      type: 'quantitative',
      axis: {
        title: 'Average of gene'
      }
    }
  }
};

// Embed the visualization in the container with id `vis`
// vegaEmbed('#vis', vlSpec);

var global_vlSpec = null

var plot_time = function(selector, data) {
  // Assign the specification to a local variable vlSpec.
  var vlSpec = {
    //"$schema": 'https://vega.github.io/schema/vega-lite/v3.json',
  $schema: 'https://vega.github.io/schema/vega-lite/v3.0.0-rc14.json',
    "title": {
      "text": data.hgnc_symbol,
      "anchor": "middle",
      "fontStyle": "italic",
      "fontWeight": "normal",
      "fontSize": 12
    },
    "width": 80,
    "height": 60,
    "data": {
      "values": data.values
    },
    "mark": 'line',
    "transform": [
      {
        "calculate": "log(datum.gene + 1) / log(2)",
        "as": "log2gene"
      }
    ],
    "encoding": {
      "x": {
        "field": "time", "type": "quantitative",
        "axis": {
          "title": "Time (h)",
          "values": [0, 24],
          "titleFont": "Helvetica Neue",
          "titleFontWeight": "normal"
        }
      },
      "y": {
        "aggregate": "average",
        "field": "log2gene",
        "type": "quantitative",
        "scale": {
          "domain": "unaggregated",
          "nice": 3,
          "zero": false
        },
        "axis": {
          "title": "Log2 TPM",
          "titleFont": "Helvetica Neue",
          "titleFontWeight": "normal"
        }
      },
      "color": {
        "title": "",
        "field": "stimulation",
        "type": "nominal",
        "scale": {
          "range": ["#FEB24C", "#E31A1C", "#800026"]
        }
      }
    }
  };

  vegaEmbed(selector, vlSpec);
}

var stat_columns = {
  't2':    'TNF at 2h vs 0h',
  't4':    'TNF at 4h vs 0h',
  't6':    'TNF at 6h vs 0h',
  't8':    'TNF at 8h  vs 0h',
  't10':   'TNF at 10h vs 0h',
  't12':   'TNF at 12h vs 0h',
  't18':   'TNF at 18h vs 0h',
  't24':   'TNF at 24h vs 0h',
  'd1':    'TNF and IL-17A (1) vs TNF',
  'd10':   'TNF and IL-17A (10) vs TNF',
  'CUX1':  'si-CUX1',
  'ELF3':  'si-ELF3',
  'LIFR':  'si-LIFR',
  'STAT3': 'si-STAT3',
  'STAT4': 'si-STAT4'
};

var scatter1 = document.getElementById('scatter1');
var state = {
  scatter1_selected: []
};
var x = null;
var xx = null;
var gene = null;

var global_scatter1_spec = null;

var table1 = null;

var main = function() {

  var scatter1_spec = {
    "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
    "data": { "values": state.stats_table },
    "selection": {
      "brush": {
        "type": "interval"
      }
    },
    "mark": "circle",
    // "transform": [
    //   {
    //     "calculate": "log(datum.t6_fold_change) / log(2)",
    //     "as": "t6_log2fc"
    //   },
    //   {
    //     "calculate": "log(datum.d10_fold_change) / log(2)",
    //     "as": "d10_log2fc"
    //   }
    // ],
    "encoding": {
      "x": {"field": "t6_log2_fold_change", "type": "quantitative"},
      "y": {"field": "d10_log2_fold_change", "type": "quantitative"},
      "color": {
        "condition": {"selection": "brush"},
        "value": "grey"
      }
    }
  };
  global_scatter1_spec = scatter1_spec;
  // vegaEmbed("#scatter1", scatter1_spec);

  // Does not work.
//  var vg_spec = vl.compile(scatter1_spec).spec;
//  var view = new vega.View(vega.parse(vg_spec))
//    .initialize(document.querySelector('#scatter1'))
//    .renderer('canvas')
//    .run();

  document.getElementById('loader').style.display = "none";

  make_plotly(state);
  
  make_plotly_buttons(state);

  make_table();
}

var make_plotly_buttons = function(state) {

  var cols = Object.keys(state.stats_table[0])
    .filter(d => d.indexOf("fold_change") != -1)

  var d = document.getElementById("scatter1-buttons");

  var r = document.createElement("div");
  r.className = "fl w-50";
  d.appendChild(r);

  var x = document.createElement("div");
  x.className = "fl pa2";
  x.innerHTML = "x-axis: ";
  r.appendChild(x);
  x = document.createElement("div");
  x.className = "fl pa2";
  var s = document.createElement("select");
  for (var i = 0; i < cols.length; i++) {
    var b = document.createElement("option");
    var label = cols[i].replace("_log2_fold_change", "");
    b.className = "fl w-20 br3 ma1";
    b.value = cols[i];
    if (label == "t6") {
      b.selected = 'selected';
    }
    b.innerHTML = label;
    s.appendChild(b);
  }
  s.onchange = function() {
    var y = this.value;
    state.scatter1_x = y;
    var rows = state.stats_table;
    // var newx = unpack(rows, state.scatter1_x).map(d => Math.log2(d));
    var newx = unpack(rows, state.scatter1_x);
    scatter1.data[0].x = newx;
    scatter1.layout.xaxis.title = stat_columns[y.replace("_log2_fold_change", "")];
    Plotly.redraw(scatter1);
  };
  x.appendChild(s);
  r.appendChild(x);

  var r = document.createElement("div");
  r.className = "fl w-50";
  d.appendChild(r);

  var x = document.createElement("div");
  x.className = "fl pa2";
  x.innerHTML = "y-axis: ";
  r.appendChild(x);
  x = document.createElement("div");
  x.className = "fl pa2";
  var s = document.createElement("select");
  for (var i = 0; i < cols.length; i++) {
    var b = document.createElement("option");
    var label = cols[i].replace("_log2_fold_change", "");
    b.className = "fl w-20 br3 ma1";
    b.value = cols[i];
    if (label == 'd10') {
      b.selected = 'selected';
    }
    b.innerHTML = label;
    s.appendChild(b);
  }
  s.onchange = function() {
    var y = this.value;
    //this.className = "br3 ma1 ba b--red";
    state.scatter1_y = y;
    var rows = state.stats_table;
    // var newy = unpack(rows, state.scatter1_y).map(d => Math.log2(d));
    var newy = unpack(rows, state.scatter1_y);
    scatter1.data[0].y = newy;
    scatter1.layout.yaxis.title = stat_columns[y.replace("_log2_fold_change", "")];
    Plotly.redraw(scatter1);
  };
  x.appendChild(s);
  r.appendChild(x);
}

function unpack(rows, key) {
  return rows.map(function(row) { return row[key]; });
}

var make_plotly = function(state) {

  var bonf = 0.05 / state.stats_table.length;

  // var rows = state.stats_table.filter(
  //   row => row.t6_pvalue <= bonf || row.d10_pvalue <= bonf
  // );
  var rows = state.stats_table;

  var data = [{
    type: 'scattergl',
    mode: 'markers',
    // x: unpack(rows, "t6_fold_change").map(d => Math.log2(d)),
    // y: unpack(rows, "d10_fold_change").map(d => Math.log2(d)),
    x: unpack(rows, "t6_log2_fold_change"),
    y: unpack(rows, "d10_log2_fold_change"),
    text: unpack(rows, 'hgnc_symbol').map(d => '<i>' + d + '</i>'),
    customdata: unpack(rows, "ensembl_id"),
    hoverinfo: 'text',
    showlegend: false,
    marker: {
      color: '#333',
      size: 6,
    },
    selected: {
      marker: {
        color: 'red'
      }
    },
    unselected: {
      marker: {
        color: '#333'
      }
    }
  }];

  var layout = {
    title: 'Log2 Fold Change',
    // dragmode: 'lasso',
    dragmode: 'select',
    hovermode: 'closest',
    hoverlabel: {
      font: {
        family: 'Helvetica Neue, Helvetica, Tahoma, sans'
      }
    },
    // showlegend: false,
    // showlegend: true,
    font: {
      size: 14,
      family: 'Helvetica Neue, Helvetica, Tahoma, sans'
    },
    xaxis: {
      title: {
        text: stat_columns["t6"]
      },
      zerolinecolor: '#969696',
      showline: true,
      linecolor: '#969696',
      linewidth: 1,
      mirror: 'ticks'
    },
    yaxis: {
      title: {
        text: 'TNF and IL-17\nvs TNF'
      },
      // tickformat: '%',
      zerolinecolor: '#969696',
      showline: true,
      linecolor: '#969696',
      linewidth: 1,
      mirror: 'ticks'
    }
  };

  var config = {
    showSendToCloud: false,
    scrollZoom: false,
    staticPlot: false,
    displaylogo: false,
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToRemove: [
      'toggleSpikelines',
      'hoverCompareCartesian', 'hoverClosestCartesian',
      'zoomIn2d', 'zoomOut2d', 'resetScale2d'
    ]
  };

  Plotly.newPlot('scatter1', data, layout, config);

  scatter1.on('plotly_hover', function(eventData) {
    if (eventData) {
      state.scatter1_hovered = eventData.points.map(
        d => rows[d.pointIndex].ensembl_id
      );
      table1.setFilter("ensembl_id", "in", state.scatter1_hovered);
    }
  });

  scatter1.on('plotly_unhover', function(eventData) {
    table1.clearFilter();
    if (state.scatter1_selected.length > 0) {
      table1.setFilter("ensembl_id", "in", state.scatter1_selected);
    }
  });

  scatter1.on('plotly_click', function(eventData) {
    // console.log(eventData.points[0].text);
  });

  scatter1.on('plotly_selected', function(eventData) {
    if (eventData) {
      console.log(eventData);
      state.scatter1_selected = eventData.points.map(
        d => rows[d.pointIndex].ensembl_id
      );
      // table1.replaceData(state.stats_table.filter(
      //    d => state.scatter1_selected.indexOf(d.ensembl_id) != -1
      // ));
      table1.setFilter("ensembl_id", "in", state.scatter1_selected);
      // var ix = eventData.points.map(d => d.pointIndex)
      // eventData.points.forEach(function(pt) {
      //   console.log(state.stats_table[pt.pointIndex]);
      // });
    }
  });

  scatter1.on('plotly_deselect', function() {
    state.scatter1_selected = [];
    table1.clearFilter();
    // table1.replaceData(state.stats_table);
    var update = {
      marker: {
        color: '#333'
      }
    }
    var trace = 0;
    // Plotly.restyle(scatter1, 'marker.color', '#333')
    Plotly.restyle(scatter1, update, trace);
  });

  // We can use this to add a new trace.
  //
  // var new_data = [{
  //   type: 'scattergl',
  //   mode: 'markers',
  //   x: [6],
  //   y: [3],
  //   text: ['butt'],
  //   customdata: ['cool'],
  //   hoverinfo: 'text',
  //   showlegend: true,
  //   name: 'butt',
  //   marker: {
  //     color: 'red',
  //     size: 6,
  //   }
  // }];
  // Plotly.addTraces(scatter1, new_data);
}



// Tabulator
// ------------------------------------------------------------------

var lineFormatter = function(cell, formatterParams, onRendered) {
  onRendered(function() {
    var el = cell.getElement();

    var div = document.createElement('div');
    div.width = 200;
    div.height = 120;
    el.appendChild(div)

    var this_ensembl_id = cell.getRow().getData().ensembl_id
    div.id = this_ensembl_id;

    var this_hgnc_symbol = cell.getRow().getData().hgnc_symbol
    gene = state.tpm_matrix.filter((d) => { return d.ID_REF == this_ensembl_id })[0];

    var m = state.metadata;
    for (var i = 0; i < m.length; i++) {
      m[i].gene = +gene[m[i].sample];
    }

    var stimulations = [
      ...new Set(m.filter(d => d.time == 2).map(d => d.stimulation))
    ];
    var blocks = m.filter(d => d.stimulation == "None");
    var new_blocks = [];
    for (var s in stimulations) {
      for (var i = 0; i < blocks.length; i++) {
        var x = Object.assign({}, blocks[i]);
        x.stimulation = stimulations[s];
        new_blocks.push(x);
      }
    }
    
    plot_time(
      '#' + this_ensembl_id,
      {
        'hgnc_symbol': this_hgnc_symbol,
        'values': m.filter(d => d.stimulation != "None").concat(new_blocks)
      }
    );
  });
};

var lineFormatter2 = function(cell, formatterParams, onRendered) {
  onRendered(function() {
    var el = cell.getElement();

    var div = document.createElement('div');
    div.width = 200;
    div.height = 120;
    el.appendChild(div)

    var this_ensembl_id = cell.getRow().getData().ensembl_id
    var this_hgnc_symbol = cell.getRow().getData().hgnc_symbol
    div.id = 'sirna-' + this_ensembl_id;

    plot_sirna(
      '#sirna-' + this_ensembl_id,
      {
        'hgnc_symbol': this_hgnc_symbol,
        'values': state.stats_table.filter(d => d.ensembl_id == this_ensembl_id)[0]
      }
    );
  });
};

var table1_columns = [
  //{title:"Ensembl ID", field:"ensembl_id", width:150},
  {
    title: "Time Series", field: "ensembl_id", width: 250, formatter: lineFormatter
  },
  {
    title: "siRNA", field: "ensembl_id", width: 150, formatter: lineFormatter2
  },
  {
    title:"Gene",
    field:"hgnc_symbol",
    headerFilter: "input",
//         headerFilterParams: {min: 0, max: 10, step: 1}
    width: 100
  },
  // {
  //   title:"Mean",
  //   field:"mean_expression",
// //         headerFilter: "numeric",
// //         headerFilterParams: {min: 0, max: 10, step: 1}
  // },
  {
    title: "TNF at 6h",
    columns: [
      {title:"log2FC", field:"t6_log2_fold_change"},
      {title:"P", field:"t6_pvalue", headerFilter: "input", headerFilterFunc: "<="}
    ]
  },
  {
    title: "IL-17A (10)",
    columns: [
      {title:"log2FC", field:"d10_log2_fold_change"},
      {title:"P", field:"d10_pvalue", headerFilter: "input", headerFilterFunc: "<="}
    ]
  },
  {
    title: "si-CUX1",
    columns: [
      {title:"log2FC", field:"CUX1_log2_fold_change"},
      {title:"P", field:"CUX1_pvalue"}
    ]
  },
  {
    title: "si-LIFR",
    columns: [
      {title:"log2FC", field:"LIFR_log2_fold_change"},
      {title:"P", field:"LIFR_pvalue"}
    ]
  },
  {
    title: "si-STAT3",
    columns: [
      {title:"log2FC", field:"STAT3_log2_fold_change"},
      {title:"P", field:"STAT3_pvalue"}
    ]
  },
  {
    title: "si-STAT4",
    columns: [
      {title:"log2FC", field:"STAT4_log2_fold_change"},
      {title:"P", field:"STAT4_pvalue"}
    ]
  },
  {
    title: "si-ELF3",
    columns: [
      {title:"log2FC", field:"ELF3_log2_fold_change"},
      {title:"P", field:"ELF3_pvalue"}
    ]
  }
];

var row_hover = function(e, row) {
  // e - the event object
  // row - row component
  var ensembl_id = row._row.data.ensembl_id;
  var point = scatter1.data[0].customdata.indexOf(ensembl_id);
  // console.log(ensembl_id);
  // console.log(point);
  console.log(scatter1.data[0]);

  // This only works for 'scatter', but not for 'scattergl'
  //
  // Plotly.Fx.hover('scatter1', [
  //   {curveNumber: 0, pointNumber: point}
  // ]);

}

var make_table = function() {
  //create Tabulator on DOM element with id "example-table"
  table1 = new Tabulator("#example-table", {
    selectable: false,
    // set height of table (in CSS or here), this enables the Virtual DOM and
    // improves render speed dramatically (can be any valid css height value)
    height: 600, 
    data: state.stats_table,
    layout: "fitDataFill",
    columns: table1_columns
    // rowMouseEnter: row_hover
  });
}

var si = null;
var plot_sirna  = function(selector, data) {
  si = data;

  var sirnas = ["CUX1", "LIFR", "STAT3", "STAT4", "ELF3"];
  var values = [];

  for (var s in sirnas) {
    values.push({
      "sirna":       sirnas[s],
      "fold_change": Math.pow(2, data.values[sirnas[s] + "_log2_fold_change"]),
      "ci_low":      Math.pow(2, data.values[sirnas[s] + "_ci_low"]),
      "ci_high":     Math.pow(2, data.values[sirnas[s] + "_ci_high"]),
      "fdr":     data.values[sirnas[s] + "_fdr"]
    })
  }

  var vlSpec = {
    //"$schema": 'https://vega.github.io/schema/vega-lite/v3.json',
  $schema: 'https://vega.github.io/schema/vega-lite/v3.0.0-rc14.json',
    "title": {
      "text": data.hgnc_symbol,
      "anchor": "middle",
      "fontStyle": "italic",
      "fontWeight": "normal",
      "fontSize": 12
    },
    "width": 80,
    "height": 60,
    "data": {
      "values": values
    },
    "transform": [
      {
        "calculate": "datum.fdr < 0.05",
        "as": "signif"
      },
      {
        "calculate": "'si-' + datum.sirna",
        "as": "sirna"
      }
    ],
    "layer": [
      {
        "data": {
          "values": [
            { "x": 1, "y": null},
            { "x": 0.5, "y": null},
            { "x": 2, "y": null}
          ]
        },
        "mark": {"type": "rect"},
        "encoding": {
          "x": {
            "field": "x", "type": "quantitative",
            "aggregate": "min",
            "axis": {
              "title": "Fold Change",
              "titleFont": "Helvetica Neue",
              "titleFontWeight": "normal"
            }
          },
          "x2": {
            "field": "x", "type": "quantitative",
            "aggregate": "max"
          },
          //"size": {"value": 1},
          "color": {"value": "#EFEFEF"}
        }
      },
      {
        "data": {
          "values": [
            { "x": 1, "y": null}
          ]
        },
        "mark": {"type": "rule"},
        "encoding": {
          "x": {
            "field": "x", "type": "quantitative"
          },
          "size": {"value": 1},
          "color": {"value": "#333"}
        }
      },
      {
        "mark": {"type": "rule"},
        "encoding": {
          "x": {
            "field": "ci_low",
            "type": "quantitative",
            "scale": {
              "zero": false
            },
            "axis": {
              "title": false
            }
          },
          "x2": {
            "field": "ci_high",
            "type": "quantitative",
            "scale": {
              "zero": false
            },
            "axis": {
              "title": false
            }
          },
          "y": {
            "field": "sirna",
            "type": "ordinal",
            "axis": {
              "title": false
            }
          },
          "color": {
            "field": "signif",
            "type": "nominal",
            "scale": {
              "range": ["#999", "#333"]
            }
          }
        }
      },
      {
        "mark": {"type": "point", "filled": true},
        "encoding": {
          "x": {
            "field": "fold_change",
            "type": "quantitative",
            "scale": {
              "zero": false
            },
            "axis": {
              "title": false
            }
          },
          "y": {
            "field": "sirna",
            "type": "ordinal",
            "axis": {
              "title": false
            }
          },
          "color": {
            "field": "signif",
            "type": "nominal",
            "scale": {
              "range": ["#999", "#333"]
            }
          }
        }
      }
    ]
  };

  vegaEmbed(selector, vlSpec);
}

fetch('data/rnaseq-data.tsv.gz', function(tsv_string) {
  state['stats_table'] = d3.tsvParse(tsv_string, function(d) {
    return {
      ensembl_id:        d.ensembl_id,
      hgnc_symbol:       d.hgnc_symbol,
      mean_expression:   +d.mean_expression,
      t2_log2_fold_change:    +d.t2_log2_fold_change,
      t2_pvalue:         +d.t2_pvalue,
      t4_log2_fold_change:    +d.t4_log2_fold_change,
      t4_pvalue:         +d.t4_pvalue,
      t6_log2_fold_change:    +d.t6_log2_fold_change,
      t6_pvalue:         +d.t6_pvalue,
      t8_log2_fold_change:    +d.t8_log2_fold_change,
      t8_pvalue:         +d.t8_pvalue,
      t12_log2_fold_change:   +d.t12_log2_fold_change,
      t12_pvalue:        +d.t12_pvalue,
      t18_log2_fold_change:   +d.t18_log2_fold_change,
      t18_pvalue:        +d.t18_pvalue,
      t24_log2_fold_change:   +d.t24_log2_fold_change,
      t24_pvalue:        +d.t24_pvalue,
      d1_log2_fold_change:    +d.d1_log2_fold_change,
      d1_pvalue:         +d.d1_pvalue,
      d10_log2_fold_change:   +d.d10_log2_fold_change,
      d10_pvalue:        +d.d10_pvalue,
      CUX1_log2_fold_change:  +d.CUX1_log2_fold_change,
      CUX1_ci_low:       +d.CUX1_ci_low,
      CUX1_ci_high:      +d.CUX1_ci_high,
      CUX1_pvalue:       +d.CUX1_pvalue,
      CUX1_fdr:          +d.CUX1_fdr,
      ELF3_log2_fold_change:  +d.ELF3_log2_fold_change,
      ELF3_ci_low:       +d.ELF3_ci_low,
      ELF3_ci_high:      +d.ELF3_ci_high,
      ELF3_pvalue:       +d.ELF3_pvalue,
      ELF3_fdr:          +d.ELF3_fdr,
      LIFR_log2_fold_change:  +d.LIFR_log2_fold_change,
      LIFR_ci_low:       +d.LIFR_ci_low,
      LIFR_ci_high:      +d.LIFR_ci_high,
      LIFR_pvalue:       +d.LIFR_pvalue,
      LIFR_fdr:          +d.LIFR_fdr,
      STAT3_log2_fold_change: +d.STAT3_log2_fold_change,
      STAT3_ci_low:      +d.STAT3_ci_low,
      STAT3_ci_high:     +d.STAT3_ci_high,
      STAT3_pvalue:      +d.STAT3_pvalue,
      STAT3_fdr:         +d.STAT3_fdr,
      STAT4_log2_fold_change: +d.STAT4_log2_fold_change,
      STAT4_ci_low:      +d.STAT4_ci_low,
      STAT4_ci_high:     +d.STAT4_ci_high,
      STAT4_pvalue:      +d.STAT4_pvalue,
      STAT4_fdr:         +d.STAT4_fdr
    };
  }).filter((d) => { return d.mean_expression > 1; });

  fetch('data/rnaseq-data-1-gene-tpm.tsv.gz', function(s) {
    state['tpm_matrix'] = d3.tsvParse(s);

    fetch('data/rnaseq-data-1-metadata.tsv.gz', function(s) {
      state['metadata'] = d3.tsvParse(s).map(d => {
        d.time = +d.time;
        return d;
      });

      main();

    })
  })
});

document.getElementById("download-select").onchange = function() {
  if (this.value == "csv") {
    table1.download("csv", "data.csv");
  }
  else if (this.value == "json") {
    table1.download("json", "data.json");
  }
  else if (this.value == "xlsx") {
    table1.download("xlsx", "data.xlsx", {sheetName:"My Data"});
  }
//  else if (this.value == "pdf") {
//    table1.download("pdf", "data.pdf", {
//      orientation: "portrait", //set page orientation to portrait
//      title: "Gene expression response to TNF and IL-17A", //add title to report
//    });
//  }
};


//define some sample data
//var tabledata = [
//   {id:1, name:"Oli Bob", age:"12", col:"red", dob:""},
//   {id:2, name:"Mary May", age:"1", col:"blue", dob:"14/05/1982"},
//   {id:3, name:"Christine Lobowski", age:"42", col:"green", dob:"22/05/1982"},
//   {id:4, name:"Brendon Philips", age:"125", col:"orange", dob:"01/08/1980"},
//   {id:5, name:"Margret Marmajuke", age:"16", col:"yellow", dob:"31/01/1999"},
//];

////create Tabulator on DOM element with id "example-table"
//var table = new Tabulator("#example-table", {
//   height:205, // set height of table (in CSS or here), this enables the Virtual DOM and improves render speed dramatically (can be any valid css height value)
//   data:tabledata, //assign data to table
//   layout:"fitColumns", //fit columns to width of table (optional)
//   columns:[ //Define Table Columns
//     {title:"Ensembl ID", field:"ensembl_id", width:150},
//     {title:"HGNC Symbol", field:"hgnc_symbol", width:150},
//     {title:"Mean", field:"mean_expression", width:150}
////     {title:"Age", field:"age", align:"left", formatter:"progress"},
////     {title:"Favourite Color", field:"col"},
////     {title:"Date Of Birth", field:"dob", sorter:"date", align:"center"},
//   ],
//   rowClick:function(e, row){ //trigger an alert message when the row is clicked
//     alert("Row " + row.getData().id + " Clicked!!!!");
//   },
//});
