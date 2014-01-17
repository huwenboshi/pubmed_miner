// check if array is already sorted
function is_table_sorted(table) {
    var sorted = true;
    for (var i = 0; i < table.length - 1; i++) {
        if (table[i].sort_key < table[i+1].sort_key) {
            sorted = false;
            break;
        }
    }
    return sorted;
}

// comparison function for sort
function compare(a,b) {
    if (a.sort_key > b.sort_key)
        return -1;
    if (a.sort_key < b.sort_key)
        return 1;
    return 0;
}

// sort a table
function sorttable(table_id) {

    // find which checkboxes are checked
    var checkboxes_status = new Array();
    var checkboxes = $("#"+table_id+" input[type='checkbox']");
    var num_checkboxes = checkboxes.length;
    var num_checked = 0;
    for(var i = 0; i < num_checkboxes; i++) {
        checkboxes_status[i] = checkboxes[i].checked;
        if (checkboxes_status[i] == true) {
            num_checked += 1;
        }
    }
    
    // iterate through table and create table content array for sorting
    var table_content = new Array();
    var table = document.getElementById(table_id);
    var num_rows = table.rows.length;
    var num_cols = table.rows[0].cells.length;
    
    for (var i = 1; i < num_rows; i++) {
        row_content = new Object();
        row_content["sort_key"] = 0;
        // skip index column
        for (var j = 1; j < num_cols; j++) {
            if (j == 1) {
                row_content["title_abstract"] = table.rows[i].cells[j].innerHTML;
            }
            else {
                row_content[j-2] = table.rows[i].cells[j].innerHTML;
                // add up sort key
                if (checkboxes_status[j-2] == true) {
                    var term_count_str = table.rows[i].cells[j].getElementsByClassName("term_count")[0].innerHTML;
                    var term_count = parseInt(term_count_str);
                    row_content["sort_key"] += term_count;
                }
            }
        }
        table_content[i-1] = row_content;
    }
    
    // apply the sort
    if (num_checked > 0 && table_content.length > 0) {
        if (!is_table_sorted(table_content)) {
            table_content.sort(compare);
        }
    }
    
    // rerender the table, skip header row
    for (var i = 1; i < num_rows; i++) {
        // skip index column
        for (var j = 1; j < num_cols; j++) {
            if (j == 1) {
                table.rows[i].cells[j].innerHTML = table_content[i-1].title_abstract;
            }
            else {
                table.rows[i].cells[j].innerHTML = table_content[i-1][j-2];
                if (checkboxes_status[j-2] == true) {
                    table.rows[i].cells[j].style.fontWeight = "bold";
                }
                else {
                    table.rows[i].cells[j].style.fontWeight = "normal";
                }
            }
        }
    }
}

// sort a table
function sortoverview(table_id) {
    // find which checkboxes are checked
    var checkboxes_status = new Array();
    var checkboxes = $("#"+table_id+" input[type='checkbox']");
    var num_checkboxes = checkboxes.length;
    var num_checked = 0;
    for(var i = 0; i < num_checkboxes; i++) {
        checkboxes_status[i] = checkboxes[i].checked;
        if (checkboxes_status[i] == true) {
            num_checked += 1;
        }
    }

    // iterate through table and prepare for sorting
    var table_content = new Array();
    var table = document.getElementById(table_id);
    var num_rows = table.rows.length;
    var num_cols = table.rows[0].cells.length;
    
    for (var i = 1; i < num_rows; i++) {
        row_content = new Object();
        row_content["sort_key"] = 0;
        // skip index column
        for (var j = 0; j < num_cols; j++) {
            // add up sort key
            if (j == 0) {
                row_content["gene_id"] = table.rows[i].cells[j].innerHTML;
            }
            else {
                row_content[j-1] = table.rows[i].cells[j].innerHTML;
                if (checkboxes_status[j-1] == true) {
                    var term_count_str = table.rows[i].cells[j].getElementsByClassName("abstract_count")[0].innerHTML;
                    var term_count = parseInt(term_count_str);
                    row_content["sort_key"] += term_count;
                }
            }
        }
        table_content[i-1] = row_content;
    }
    
    // apply the sort
    if (num_checked > 0 && table_content.length > 0) {
        if (!is_table_sorted(table_content)) {
            table_content.sort(compare);
        }
    }
    
    // rerender the table, skip header row
    for (var i = 1; i < num_rows; i++) {
        // skip index column
        for (var j = 0; j < num_cols; j++) {
            if (j == 0) {
                table.rows[i].cells[j].innerHTML = table_content[i-1]["gene_id"];
            }
            else {
                table.rows[i].cells[j].innerHTML = table_content[i-1][j-1];
                if (checkboxes_status[j-1] == true) {
                    table.rows[i].cells[j].style.fontWeight = "bold";
                }
                else {
                    table.rows[i].cells[j].style.fontWeight = "normal";
                }
            }
        }
    }
}
