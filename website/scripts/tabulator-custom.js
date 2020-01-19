calc_significant = function(values, data, calcParams) {
  //values - array of column values
  //data - all table data
  //calcParams - params passed from the column definition object

  var count = 0;

  values.forEach(function(value) {
    if (value < 0.05) {
      count++;
    }
  });

  return count;
}
