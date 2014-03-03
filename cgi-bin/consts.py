# tool, email info
tool_email = 'tool=pubmed_article_rec&email=pubmed_article_rec@yahoo.com'

# params
dbfrom_db_lknm = '&dbfrom=gene&db=pubmed&linkname=gene_pubmed&'
dbfrom_db_lknm_cmd = dbfrom_db_lknm+'cmd=neighbor_history'
efetch_fixed_param = '&db=pubmed&rettype=abstract&retmode=xml'
query_param = 'db=gene&cmd=Retrieve&dopt=full_report&list_uids='

# urls 
pubmed_url = 'http://www.ncbi.nlm.nih.gov/pubmed'
eutils_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils'
esummary_url = eutils_url+'/esummary.fcgi?'+tool_email
elink_url = eutils_url+'/elink.fcgi?'+tool_email
efetch_url = eutils_url+'/efetch.fcgi?'+tool_email
gene_web_url = 'http://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&dopt='
gene_web_url += 'full_report&list_uids='
genomernai_url = 'http://genomernai.de/GenomeRNAi/genedetails/'
biogpsorg_url = 'http://biogps.org/#goto=genereport&id='
genecards_url = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene='
ucsc_url = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&position='

# http header
http_header = 'Content-Type: text/html Connection: keep-alive'

# html header
html_header = """
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" 
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">
    <head>
    <meta charset="utf-8"/>
    <title>Search Result</title>
    <script src="http://code.jquery.com/jquery-1.10.1.min.js"></script>
    <script src="http://code.jquery.com/jquery-migrate-1.2.1.min.js"></script>
    <script src="../javascripts/sorttable.js"></script>
    <script type="text/javascript">
        $(document).ready(function(){
            $('.local_opt').live('change', function() {
                var id_str = $(this).attr('id');
                var id_arr = id_str.split('_');
                var gene_id = id_arr[4];
                var table_id = "articles_gene_id_"+gene_id;
                sorttable(table_id);
            });
            $('.overview_opt').live('change', function() {
                sortoverview("overview_top");
            });
            $('.show_more_less').live('click', function() {
                $div = $(this).parent();
                if ($div.data('open')) {
                    $div.css({height:'75px',overflow:'hidden',
                        'line-height':'25px'});
                    $div.data('open', 0);
                    $(this).html("show more");
                }
                else {
                    $div.css({height:'100%'});
                    $div.data('open', 1);
                    $(this).html("show less");
                }
            });
            $('.show_hide').live('click', function() {
                if ($(this).html() == "hide") {
                    $(this).siblings("table").css({display:'none'});
                    $(this).html("show");
                }
                else {
                    $(this).siblings("table").css({display:'inline-table'});
                    $(this).html("hide");
                }
            });
            $("#overview_top").html($("#overview_bottom").html());
            $("#overview_bottom").html("");
            $("#loading").css({display: 'none'});
            
            $("#gene_term_network_top").
                html($("#gene_term_network_bottom").html());
            $("#gene_term_network_bottom").html("");
            $("#loading_gene_term_network").css({display: 'none'});
            
            // Function to get the max value in an Array
            Array.max = function(array){
                return Math.max.apply(Math,array);
            };
         
            // Get all data values from our table cells making
            // sure to ignore the first column of text
            // Use the parseInt function to convert the text string to a number
         
            var counts= $('.heat-map tbody td').not('.stats-title').map(function() {
                return parseInt($(this).text());
            }).get();
            
            // run max value function and store in variable
            var max = Array.max(counts);
            n = 100; // Declare the number of groups
         
            // Define the ending colour, which is white
            xr = 255; // Red value
            xg = 255; // Green value
            xb = 255; // Blue value
         
            // Define the starting colour #f32075
            yr = 52; // Red value
            yg = 119; // Green value
            yb = 220; // Blue value
            
            // Loop through each data point and calculate its % value
            $('.heat-map tbody td').not('.stats-title').each(function(){
                var val = parseInt($(this).text());
                var pos = parseInt((Math.round((val/max)*100)).toFixed(0));
                red = parseInt((xr + (( pos * (yr - xr)) / (n-1))).toFixed(0));
                green = parseInt((xg + (( pos * (yg - xg)) / (n-1))).toFixed(0));
                blue = parseInt((xb + (( pos * (yb - xb)) / (n-1))).toFixed(0));
                clr = 'rgb('+red+','+green+','+blue+')';
                $(this).css({backgroundColor:clr});
            });
        });
    </script>
    <style>
        table {
            border-bottom: 1px Solid Black;         
            border-right: 1px Solid Black;         
            border-collapse : collapse;  
        }
        table td, table th {    
            border-left: 1px Solid Black;         
            border-top: 1px Solid Black;              
            border-bottom: 1px Solid Black;    
            border-right:none;  
        }
        table td {
            max-width:800px;
            word-wrap:break-word;
        }
        .gwas_title {
            max-width:450px;
            word-wrap:break-word;
        }
        .gwas_trait {
            max-width:150px;
            word-wrap:break-word;
        }
        .gwas_snp {
            max-width:100px;
            word-wrap:break-word;
        }
        .gwas_gene {
            max-width:150px;
            word-wrap:break-word;
        }
        .abstract_txt {
            line-height:25px;
            height:75px;
            overflow: hidden;
         }
        #overview_bottom {
            display: none;
        }
    </style>
    </head>
"""
