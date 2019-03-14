
var fetch = function(data_url, callback) {
  logTime('downloading data');
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
      "values": data.values.filter((d) => {return d.gene > 0})
    },
    "mark": 'line',
    "transform": [
      {"calculate": "log(datum.gene + 1) / log(2)", "as": "log2gene"}
    ],
    "encoding": {
      "x": {
        "field": "time", "type": "quantitative",
        "axis": {
          title: false
        }
      }, //"type": "ordinal", "timeUnit": "hours"},
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
        "field": "stimulation",
        "type": "nominal"
        // "type": "ordinal"
      }
    }
  };
  global_vlSpec = vlSpec;

  // Embed the visualization in the container with id `vis`
  // vegaEmbed('#vis', vlSpec);
  vegaEmbed(selector, vlSpec);
}

var state = {};
var x = null;
var gene = null;

var main = function() {

  var lineFormatter = function(cell, formatterParams, onRendered){
    onRendered(function() { //instantiate sparkline after the cell element has been aded to the DOM
      el = cell.getElement();

      //var canvas = document.createElement('canvas');
      //canvas.id = "CursorLayer";
      //canvas.width = 200;
      //canvas.height = 100;
      //canvas.style.zIndex = 8;
      //canvas.style.position = "absolute";
      //canvas.style.border = "1px solid";
      //el.appendChild(canvas);
      
      var div = document.createElement('div');
      div.width = 200;
      div.height = 100;
      el.appendChild(div)
      var this_ensembl_id = cell.getRow().getData().ensembl_id
      div.id = this_ensembl_id;
      var this_hgnc_symbol = cell.getRow().getData().hgnc_symbol
      gene = state.tpm_matrix.filter((d) => { return d.ID_REF == this_ensembl_id })[0];
      x = Object.assign({}, state.metadata);
      for (var i = 0; i < 175; i++) {
        x[i].gene = +gene[x[i].sample];
        x[i].time = +x[i].time;
      }
      var retval = [];
      for (var i in x) {
        if (x[i].stimulation != "None") {
          retval.push(x[i]);
        }
      }
      x = retval;
      update_chart('#' + this_ensembl_id, {'hgnc_symbol': this_hgnc_symbol, 'values': x});
    });
  };

  //create Tabulator on DOM element with id "example-table"
  var table = new Tabulator("#example-table", {
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
     columns: [ //Define Table Columns
       //{title:"Ensembl ID", field:"ensembl_id", width:150},
       {title:"Plot", field: "ensembl_id", width: 250, formatter: lineFormatter},
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
       }
     ],
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

fetch('data/rnaseq-data-1.tsv.gz', function(tsv_string) {
  state['stats_table'] = d3.tsvParse(tsv_string, function(d) {
    return {
      ensembl_id:      d.ensembl_id,
      hgnc_symbol:     d.hgnc_symbol,
      mean_expression: +d.mean_expression,
      t6_fold_change:  +d.t6_fold_change,
      t6_pvalue:       +d.t6_pvalue,
      d10_fold_change: +d.d10_fold_change,
      d10_pvalue:      +d.d10_pvalue
    };
  }).filter((d) => { return d.mean_expression > 1; });

  fetch('data/rnaseq-data-1-gene-tpm.tsv.gz', function(s) {
    state['tpm_matrix'] = d3.tsvParse(s);

    fetch('data/rnaseq-data-1-metadata.tsv.gz', function(s) {
      state['metadata'] = d3.tsvParse(s);

      main();
      // update_chart('#vis', {'values': mydata});
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
