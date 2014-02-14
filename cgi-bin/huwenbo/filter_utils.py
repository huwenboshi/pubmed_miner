# display ewas query result
def print_ewas_query_result(ewas_query_result):
    
    ########## gene table ##########
    gene_exp_result = ewas_query_result[0]
    print '<h3>Genes implicated by CG methylations '
    print 'associated with gene expression</h3>'
    print '<table>'
    print """<tr>
                <th>mCG genome-wide<br/> position</th>
                <th>1Mb-window</th>
                <th>Implicated mouse<br/>gene probe ID</th>
                <th>Gene position</th>
                <th>p-value</th>
                <th>Human ortholog<br/>entrez ID</th>
            </tr>"""
    for result in gene_exp_result:
        print '<tr>'
        
        # mCG 
        print '<td>%s</td>' % result[7]
        
        # mCG window
        window_st = result[6]-500000 if(result[6]-500000 > 0) else 0
        window_ed = result[6]+500000
        window_str = "chr%s:%s-%s" % (result[5], str(window_st),
            str(window_ed))
        print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
            window_str, window_str)
        
        # probe ID
        print '<td>%s</td>' % result[15]
        
        # gene position
        gene_pos_str = 'chr%s:%s-%s' % (result[10], result[11], result[12])
        print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
            gene_pos_str, gene_pos_str)
        
        # p-value, human entrez gene id
        print '<td>%s</td><td>%s</td>' % (result[0], result[18])
        print '</tr>'
    print '</table>'
    
    print '<br/>'
    print '<hr/>'
    
    ########## protein table ##########
    prot_exp_result = ewas_query_result[1]
    print '<h3>Genes implicated by CG methylations associated '
    print 'with protein expression</h3>'
    print '<table>'
    print """<tr>
                <th>mCG genome-wide<br/> position</th>
                <th>1Mb-window</th>
                <th>Implicated mouse<br/>gene entrez ID</th>
                <th>Gene position</th>
                <th>p-value</th>
                <th>Human ortholog<br/>entrez ID</th>
            </tr>"""
    for result in prot_exp_result:
        print '<tr>'
        # mCG 
        print '<td>%s</td>' % result[7]
        
        # mCG window
        window_st = result[6]-500000 if(result[6]-500000 > 0) else 0
        window_ed = result[6]+500000
        window_str = "chr%s:%s-%s" % (result[5], str(window_st),
            str(window_ed))
        print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
            window_str, window_str)
        
        # entrez id
        print '<td>%s</td>' % result[15]
        
        # gene position
        gene_pos_str = 'chr%s:%s-%s' % (result[10], result[11], result[12])
        print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
            gene_pos_str, gene_pos_str)
        
        # p-value, human entrez gene id
        print '<td>%s</td><td>%s</td>' % (result[0], result[18])
        print '</tr>'
    print '</table>'
    
    print '<br/>'
    print '<hr/>'
    
    ########## trait table ##########
    trait_exp_result = ewas_query_result[2]
    print '<h3>Genes implicated by CG methylations associated '
    print 'with clinical and metabolite trait</h3>'
    print '<table>'
    print """<tr>
                <th>mCG genome-wide<br/>position</th>
                <th>1-Mb window</th>
                <th>Implicated mouse<br/>gene transcript ID</th>
                <th>Gene position</th>
                <th>Phenotype</th>
                <th>Phenotype class</th>
                <th>p-value</th>
                <th>Human ortholog<br/>entrez ID</th>
            </tr>"""
    for result in trait_exp_result:
        print '<tr>'
        
        # mCG position
        print '<td>%s</td>' % result[6]
        
        # mCG window
        window_st = result[5]-500000 if(result[5]-500000 > 0) else 0
        window_ed = result[5]+500000
        window_str = "chr%s:%s-%s" % (result[4], str(window_st),
            str(window_ed))
        print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
            window_str, window_str)
        
        # transcript id
        print '<td>%s</td>' % result[17]
        
        # gene position
        gene_pos_str = 'chr%s:%s-%s' % (result[12], result[13], result[14])
        print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
            gene_pos_str, gene_pos_str)
        
        print '<td>%s</td><td>%s</td>' % (result[2], result[3])
        print '<td>%s</td><td>%s</td>' % (result[0], result[20])
        print '</tr>'
    print '</table>'
