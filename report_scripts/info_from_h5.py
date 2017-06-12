
def h5_parameter_table(h5_file, name):
    from ac_analysis.model.annotations import load_from_h5
    h5 =  load_from_h5(h5_file)

    token = str(h5.get_parameters().chip_init['exp_token'])
    experiment_mode =  h5.get_parameters().chip_init['experiment_mode']
    notes =  h5.get_parameters().chip_init['notes']
    questions_unicode =  h5.get_parameters().chip_init['question']
    questions = questions_unicode.decode('unicode_escape').encode('ascii','ignore')[1:-1]

    token_html = ('<a href = "http://10.21.53.10/tokens?view=detail&id=" %s> %s </a>' % (token,token))


    tbl = [['Token:', token_html],
           ['Mode: ', experiment_mode],
           ['Notes: ', notes],
           ['Questions:',questions]]

    ## code for table to html
    cols = ["<td>{0}</td>".format("</td><td>".join(t)) for t in tbl]

    # then use it to join the rows (tr)
    rows = "<tr>{0}</tr>".format("</tr>\n<tr>".join(cols))

    html = ("""<HTML>
                     <h1><font size ="3"> %s  </font size></h1 > 
                        <table>
                         {0}
                       </table><br>
    
                      </HTML>""".format(rows)
                        % str(name))
    return html

#
# with open(('try.html'), 'w') as output:
#     output.write(h5_parameter_table('annotations.h5'))
#
