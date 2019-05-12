
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
  $schema: 'https://vega.github.io/schema/vega-lite/v3.json',
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

var update_chart = function(selector, data) {
  // Assign the specification to a local variable vlSpec.
  var vlSpec = {
    "$schema": 'https://vega.github.io/schema/vega-lite/v3.json',
//    "title": {
//      "text": data.hgnc_symbol,
//      "anchor": "middle",
//      "fontStyle": "italic",
//      "fontSize": 16
//    },
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
          title: false,
          values: [0, 24]
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
          // "title": "Log2 TPM"
          title: false
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
  't2':    'TNF (2h)',
  't4':    'TNF (4h)',
  't6':    'TNF (6h)',
  't8':    'TNF (8h)',
  't10':   'TNF (10h)',
  't12':   'TNF (12h)',
  't18':   'TNF (18h)',
  't24':   'TNF (24h)',
  'd1':    'IL-17A (1)',
  'd10':   'IL-17A (10)',
  'CUX1':  'si-CUX1',
  'ELF3':  'si-ELF3',
  'LIFR':  'si-LIFR',
  'STAT3': 'si-STAT3',
  'STAT4': 'si-STAT4'
};

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
    "transform": [
      {
        "calculate": "log(datum.t6_fold_change) / log(2)",
        "as": "t6_log2fc"
      },
      {
        "calculate": "log(datum.d10_fold_change) / log(2)",
        "as": "d10_log2fc"
      }
    ],
    "encoding": {
      "x": {"field": "t6_log2fc", "type": "quantitative"},
      "y": {"field": "d10_log2fc", "type": "quantitative"},
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

  make_plotly(state);

  make_table();
}

function unpack(rows, key) {
  return rows.map(function(row) { return row[key]; });
}

var make_plotly = function(state) {

  var bonf = 0.05 / state.stats_table.length;

  var rows = state.stats_table.filter(
    row => row.t6_pvalue <= bonf || row.d10_pvalue <= bonf
  );

  var data = [{
    type: 'scatter',
    mode: 'markers',
    x: unpack(rows, "t6_fold_change").map(d => Math.log2(d)),
    y: unpack(rows, "d10_fold_change").map(d => Math.log2(d)),
    text: unpack(rows, 'hgnc_symbol').map(d => '<i>' + d + '</i>'),
    hoverinfo: 'text',
    marker: {
      color: '#333',
      size: 5
    },
    transforms: [
      // {
      //   type: 'filter',
      //   target: unpack(state.stats_table, 't6_pvalue'),
      //   operation: '<',
      //   value: '0.05'
      // }
    ]
  }];

  var layout = {
    // title: 'Genes',
    // dragmode: 'lasso',
    dragmode: 'select',
    hovermode: 'closest',
    hoverlabel: {
      font: {
        family: 'Helvetica Neue, Helvetica, Tahoma, sans'
      }
    },
    showlegend: false,
    font: {
      size: 14,
      family: 'Helvetica Neue, Helvetica, Tahoma, sans'
    },
    xaxis: {
      title: {
        text: 'Log2 fold-change for TNF (6h)'
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

  var scatter1 = document.getElementById('scatter1');

  Plotly.newPlot('scatter1', data, layout, config);

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
  });
}



// Tabulator
// ------------------------------------------------------------------

var lineFormatter = function(cell, formatterParams, onRendered) {
  onRendered(function() {
    var el = cell.getElement();

    var div = document.createElement('div');
    div.width = 200;
    div.height = 100;
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
    
    update_chart(
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
    div.height = 100;
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
    title: "Time Series", field: "ensembl_id", width: 240, formatter: lineFormatter
  },
  {
    title: "siRNA", field: "ensembl_id", width: 140, formatter: lineFormatter2
  },
  {
    title:"Gene",
    field:"hgnc_symbol",
    headerFilter: "input",
//         headerFilterParams: {min: 0, max: 10, step: 1}
    width: 100
  },
  {
    title:"Mean",
    field:"mean_expression",
//         headerFilter: "numeric",
//         headerFilterParams: {min: 0, max: 10, step: 1}
  },
  {
    title: "TNF",
    columns: [
      {title:"FC", field:"t6_fold_change"},
      {title:"P", field:"t6_pvalue"}
    ]
  },
  {
    title: "IL-17A",
    columns: [
      {title:"FC", field:"d10_fold_change"},
      {title:"P", field:"d10_pvalue"}
    ]
  },
  {
    title: "si-CUX1",
    columns: [
      {title:"FC", field:"CUX1_fold_change"},
      {title:"P", field:"CUX1_pvalue"}
    ]
  },
  {
    title: "si-LIFR",
    columns: [
      {title:"FC", field:"LIFR_fold_change"},
      {title:"P", field:"LIFR_pvalue"}
    ]
  },
  {
    title: "si-STAT3",
    columns: [
      {title:"FC", field:"STAT3_fold_change"},
      {title:"P", field:"STAT3_pvalue"}
    ]
  },
  {
    title: "si-STAT4",
    columns: [
      {title:"FC", field:"STAT4_fold_change"},
      {title:"P", field:"STAT4_pvalue"}
    ]
  },
  {
    title: "si-ELF3",
    columns: [
      {title:"FC", field:"ELF3_fold_change"},
      {title:"P", field:"ELF3_pvalue"}
    ]
  }
];


var make_table = function() {
  //create Tabulator on DOM element with id "example-table"
  table1 = new Tabulator("#example-table", {
    selectable: false,
     height: 800, // set height of table (in CSS or here), this enables the Virtual DOM and improves render speed dramatically (can be any valid css height value)
     data: state['stats_table'],
     // layout: "fitColumns", //fit columns to width of table (optional)
     layout: "fitDataFill",
//     rowFormatter:function(row, data){
//        //row - JQuery object for row
//        //data - the data for the row
//
//        row.css({"height":"50px"});
//    },
     columns: table1_columns,
//     rowClick: function(e, row){ //trigger an alert message when the row is clicked
//       this_ensembl_id = row.getData().ensembl_id;
//       console.log("Clicked " + this_ensembl_id);
//       gene = state.tpm_matrix.filter((d) => { return d.ID_REF == this_ensembl_id })[0];
//       x = Object.assign({}, state.metadata);
//       for (var i = 0; i < 175; i++) {
//         x[i].gene = +gene[x[i].sample];
//         x[i].time = +x[i].time;
//       }
//       var retval = [];
//       for (var i in x) {
//         if (x[i].stimulation != "None") {
//           retval.push(x[i]);
//         }
//       }
//       x = retval;
//       update_chart('#vis', {'hgnc_symbol': row.getData().hgnc_symbol, 'values': x});
//       console.log(gene);
//     }
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
      "fold_change": data.values[sirnas[s] + "_fold_change"],
      "ci_low":      data.values[sirnas[s] + "_ci_low"],
      "ci_high":     data.values[sirnas[s] + "_ci_high"],
      "fdr":     data.values[sirnas[s] + "_fdr"]
    })
  }

  var vlSpec = {
    "$schema": 'https://vega.github.io/schema/vega-lite/v3.json',
//    "title": {
//      "text": data.hgnc_symbol,
//      "anchor": "middle",
//      "fontStyle": "italic",
//      "fontSize": 16
//    },
    "width": 80,
    "height": 60,
    "data": {
      "values": values
    },
    "transform": [
      {
        "calculate": "datum.fdr < 0.05",
        "as": "signif"
      }
    ],
    "layer": [
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
      t2_fold_change:    +d.t2_fold_change,
      t2_pvalue:         +d.t2_pvalue,
      t4_fold_change:    +d.t4_fold_change,
      t4_pvalue:         +d.t4_pvalue,
      t6_fold_change:    +d.t6_fold_change,
      t6_pvalue:         +d.t6_pvalue,
      t8_fold_change:    +d.t8_fold_change,
      t8_pvalue:         +d.t8_pvalue,
      d1_fold_change:    +d.d1_fold_change,
      d1_pvalue:         +d.d1_pvalue,
      d10_fold_change:   +d.d10_fold_change,
      d10_pvalue:        +d.d10_pvalue,
      CUX1_fold_change:  +d.CUX1_fold_change,
      CUX1_ci_low:       +d.CUX1_ci_low,
      CUX1_ci_high:      +d.CUX1_ci_high,
      CUX1_pvalue:       +d.CUX1_pvalue,
      CUX1_fdr:          +d.CUX1_fdr,
      ELF3_fold_change:  +d.ELF3_fold_change,
      ELF3_ci_low:       +d.ELF3_ci_low,
      ELF3_ci_high:      +d.ELF3_ci_high,
      ELF3_pvalue:       +d.ELF3_pvalue,
      ELF3_fdr:          +d.ELF3_fdr,
      LIFR_fold_change:  +d.LIFR_fold_change,
      LIFR_ci_low:       +d.LIFR_ci_low,
      LIFR_ci_high:      +d.LIFR_ci_high,
      LIFR_pvalue:       +d.LIFR_pvalue,
      LIFR_fdr:          +d.LIFR_fdr,
      STAT3_fold_change: +d.STAT3_fold_change,
      STAT3_ci_low:      +d.STAT3_ci_low,
      STAT3_ci_high:     +d.STAT3_ci_high,
      STAT3_pvalue:      +d.STAT3_pvalue,
      STAT3_fdr:         +d.STAT3_fdr,
      STAT4_fold_change: +d.STAT4_fold_change,
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
