#!/usr/bin/env python

PCAT_LOGO = """
#######################################################################
#               _______     ______       _     _________              #
#               |_   __ \ .' ___  |     / \   |  _   _  |             #
#               | |__) |/ .'   \_|     / _ \  |_/ | | \_|             #
#               |  ___/ | |           / ___ \     | |                 #
#               _| |_    \ `.___.'\ _/ /   \ \_  _| |_                #
#              |_____|    `.____ .'|____| |____||_____|               #
#                                                                     #
#                       _,'|             _.-''``-...___..--';)        #
#                      /_ \'.      __..-' ,      ,--...--'''           #
#                    <\    .`--'''       `     /'                     #
#                      `-';'               ;   ; ;                    #
#                __...--''     ___...--_..'  .;.'                     #
#               (,__....----'''       (,..--''                        #
#                                                                     #
#######################################################################
"""

def main():
    import os
    import os.path
    from string import join
    
    BASE_PATH = os.path.expanduser("~/public_html/summaries/")
    file_list = os.listdir(BASE_PATH)
    
    LOCKED_PLOT_WIDTH = 1530//1.5
    LOCKED_PLOT_HEIGHT = 146//1.5
    
    htmls = []
    for element in file_list:
    	if (".html" in element) and element != "index.html":
    		htmls.append(element)
    
    
    html_body = """
    <html>
    <head>
    	<title>PCAT Summary Pages</title>\
    	<link rel="stylesheet" type="text/css" href="../../style/main.css">\
    </head>
    <body>
      <span>
    	<table style="text-align: left; width: 1100; height: 100   px; margin-left:auto; margin-right: auto;" border="1" cellpadding="1" cellspacing="1">
           <tbody>
            <tr>
            <td style="text-align: center; vertical-align: top; background-color:SpringGreen;"><big><big><big><big><span style="font-weight: bold;">PCAT Summary Pages</span></big></big></big></big>
            </td>
            </tr>
            </tbody>
         </table>
    	 </br><div align='center'><pre>{0}</pre></div></br>
    <table border="1" style="text-align: center; width: 1200; height: 67px; margin-left:auto; margin-right: auto; background-color: white;" border="1" cellpaddi\
    ng="2" cellspacing="2" align=center><col width=250> <col width=400>
    
    """.format(PCAT_LOGO)
    
    # Sort the html list
    htmls.sort()
    #print htmls
    
    
    #TODO: ADD GPS TIMES IN COLUMN, HEADER: <th align='center'><big><big>GPS interval</big></big></th> TABLE ITEM: <td>GPSTIMESTUFF<td>
    html_body += "<tr><th align='center'><big><big>Summary Page Link</big></big></th><th align='center'><big><big>Segment Information</big></big></th></tr>"
    for index, summarypage in enumerate(htmls):
    	summary_name = join(summarypage.split(".")[:-1])
    	plot_name = "lock-plot_{0}.png".format(summary_name)
    	#with open() as times_list:
    	
    	
    	html_body += "<tr><td><a href='./{0}'>{1}</a></td><td align='center'><img src='./img/{2}' width='{3}' height='{4}' align='center'></img></tr>" .format(summarypage, summary_name, plot_name, LOCKED_PLOT_WIDTH,LOCKED_PLOT_HEIGHT)
    
    # Close the table
    html_body += "</table>"
    
    
    # Close the html document tags
    html_body += "</span></body></html>"
    
    # Write file:
    index_html = open(os.path.expanduser("~/public_html/summaries/index.html"), "w")
    print >>index_html, html_body
    index_html.close()

# Done./end

if __name__ == '__main__':
    main()
