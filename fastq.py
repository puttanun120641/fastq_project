import gzip
import re, numpy
from Bio import SeqIO
import pandas as pd
import numpy as np

bc01 = {}
filename = 0
df = 0

def extractData(filenames):
    global filename
    filename = filenames
    extractdata()

def extractdata():
    global bc01
    global filename
    global df
    with gzip.open(filename, 'rt') as f:
        for idx, seq_record in enumerate(SeqIO.parse(f, "fastq")):
            bc01[seq_record.id] = seq_record
            # if idx == 100:
            #     break
        extractdata = []
        #writer = csv.writer(csvOut)
        for seq_id, seq_record in bc01.items():
            content = str(seq_record.description)
            runid = re.search(r'runid=(\w+)', content)
            barcode = re.search(r'barcode=(\S+)',content)
            channel = re.search(r'ch=([\w]+)', content)
            extractdata.append([seq_record.id, runid.group(1), len(seq_record.seq), numpy.median(seq_record.letter_annotations["phred_quality"]),numpy.mean(seq_record.letter_annotations["phred_quality"]) ,channel.group(1), barcode.group(1)])
    #print(extractdata)
    df = pd.DataFrame(extractdata, columns=['id', 'runid', 'Readlenght', 'Phred_quality', 'MeanPhred_quality','channel','barcode'])
    df_long = df.melt(id_vars='Readlenght',
                             var_name='id',                              
                             value_name='quality')
    df_long.head()

def calculate_N50(list_of_lengths):
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
    return median

def summaryTable(passreadscore):
    medlen = np.median(df["Readlenght"]) #medreadlenght
    medphred = np.median(df["Phred_quality"]) #medphredscore
    passread = df[(df.Phred_quality > passreadscore)] #status passread
    allread = df[(df.Phred_quality > 0)] #status allread this is df
    df["Readlenght"].sum() #total base
    passdf = passread["Readlenght"].sum()
    alldf = allread["Readlenght"].sum()
    passmedlen = np.median(passread["Readlenght"])
    passmedphred = np.median(passread["Phred_quality"])
    N50_set_allread = sorted(allread["Readlenght"])
    N50_set_passread = sorted(passread["Readlenght"])
    number_barcode_allread = np.unique(allread["barcode"])
    number_barcode_passread = np.unique(passread["barcode"])
    runid_allread = np.unique(allread["runid"])
    runid_passread = np.unique(passread["runid"])
    general_data = []
    general_data.append(["All Reads",len(runid_allread),len(number_barcode_allread)])
    general_data.append(["Pass Reads",len(runid_passread),len(number_barcode_passread)])
    general_df = pd.DataFrame(general_data, columns=['Status', 'Number of Runids', 'Number of Barcodes'])
    print(general_df)
    #basecall summary dataframe
    basecall_data = []
    basecall_data.append(["All Reads",len(allread),alldf,calculate_N50(N50_set_allread),medlen,medphred])
    basecall_data.append(["Pass Reads",len(passread),passdf,calculate_N50(N50_set_passread),passmedlen,passmedphred])
    basecall_df = pd.DataFrame(basecall_data, columns=['Status', 'Reads','Bases','N50','Median Read Length','Median PHRED Score'])
    print(basecall_df)





def figure():
    medlen = np.median(df["Readlenght"]) #medreadlenght
    medphred = np.median(df["Phred_quality"]) #medphredscore
    passread = df[(df.Phred_quality > 22)] #status passread
    allread = df[(df.Phred_quality > 0)] #status allread this is df
    df["Readlenght"].sum() #total base
    passdf = passread["Readlenght"].sum()
    alldf = allread["Readlenght"].sum()
    passmedlen = np.median(passread["Readlenght"])
    passmedphred = np.median(passread["Phred_quality"])
    N50_set_allread = sorted(allread["Readlenght"])
    N50_set_passread = sorted(passread["Readlenght"])
    number_barcode_allread = np.unique(allread["barcode"])
    number_barcode_passread = np.unique(passread["barcode"])
    runid_allread = np.unique(allread["runid"])
    runid_passread = np.unique(passread["runid"])
    #general run summary
    general_data = []
    general_data.append(["All Reads",len(runid_allread),len(number_barcode_allread)])
    general_data.append(["Pass Reads",len(runid_passread),len(number_barcode_passread)])
    general_df = pd.DataFrame(general_data, columns=['Status', 'Number of Runids', 'Number of Barcodes'])
    #print(general_df)
    #basecall summary dataframe
    basecall_data = []
    basecall_data.append(["All Reads",len(allread),alldf,calculate_N50(N50_set_allread),medlen,medphred])
    basecall_data.append(["Pass Reads",len(passread),passdf,calculate_N50(N50_set_passread),passmedlen,passmedphred])
    basecall_df = pd.DataFrame(basecall_data, columns=['Status', 'Reads','Bases','N50','Median Read Length','Median PHRED Score'])


    from pretty_html_table import build_table
    html_table1 = build_table(general_df, 'grey_dark',font_size='medium',font_family='Open Sans, sans-serif'
                            , text_align='center',even_color='black',width='auto')
    html_table2 = build_table(basecall_df, 'grey_dark',font_size='medium',font_family='Open Sans, sans-serif'
                            , text_align='center',even_color='black',width='auto')
    table1_html = '<div id="table1"><h3>General run summary</h3>'+html_table1 +'</div>'
    table1_html.replace("<table",'<table class="table"')
    table2_html = '<div id="table1"><h3>Basecall summary</h3>'+ html_table2+'</div>'
    table2_html.replace("<table",'<table class="table"')
    #filter using pandas
    passread = df[(df.Phred_quality > 22)]
    #print(passread)
    passread.size
    b = df[(df.Readlenght > 100)]
    #print(b)
    ###############################GRAPH    PLOTTING#####################################
    import seaborn as sns
    import matplotlib.pyplot as plt
    '''embended graph in to html'''
    import io
    import base64
    def fig_to_svg(fig):
        img = io.StringIO()
        fig.savefig(img, format='svg')
        return img.getvalue()
    def fig_to_png_base64(fig):
        img = io.BytesIO()
        fig.savefig(img, format='png',
                    bbox_inches='tight')
        img.seek(0)
        return base64.b64encode(img.getvalue())
    ##Kernel density estimation -density plot
    #1.Kernel density estimation- Readlenght
    '''Readlenght density'''
    sns.displot(x="Readlenght", kind="kde" ,rug=True, data =df, height=6, aspect=9/6)
    plt.xlabel("Readlenght", size=14)
    #plt.ylabel("density", size=14)
    plt.axvline(x=df.Readlenght.mean(),ls='--',color='red', label='Mean = {}'.format(df.Readlenght.mean()))
    plt.axvline(x=df.Readlenght.median(),ls='--', color='black', label='Median = {}'.format(df.Readlenght.median()))
    plt.legend(loc='upper right')
    plt.xlabel("Readlenght", size=14)
    plt.ylabel("Density", size=14)
    plt.grid()
    fig = plt.gcf()
    ############### Embended graph to html ######################
    p3_encoded_svg = fig_to_svg(fig)
    p3_html_svg = '<dir id="p2"><h3>Basecalled Readlenght</h3>{}<h3>Pass reads<h3></dir>'.format(p3_encoded_svg)
    plt.show()
    plt.clf()

    #2.Kernel density estimation- MeanPhred_quality
    '''MeanPhred_quality density'''
    sns.displot(df, x="Phred_quality", rug=True, kind="kde",height=6, aspect=9/6)
    plt.xlabel("PHRED quality scores", size=14)
    plt.ylabel("Density", size=14)
    plt.axvline(x=df.MeanPhred_quality.mean(),ls='--',color='red', label='Mean = {}'.format(df.MeanPhred_quality.mean()) )
    plt.axvline(x=medphred,ls='--',color='black', label='Median = {}'.format(medphred))
    #plt.title(' MeanPhred_quality ')
    plt.legend(loc='upper right')
    plt.grid()
    fig = plt.gcf()
    p4_encoded_svg = fig_to_svg(fig)
    p4_html_svg = '<dir id="p2"><h3>Basecalled reads PHRED quality</h3>{}<h3>Pass reads<h3></dir>'.format(p4_encoded_svg)
    plt.show()
    plt.clf()

    ## readlen&meanQ
    plt.figure()
    sns.kdeplot(x=df.Readlenght, y=df.MeanPhred_quality, cbar=True, cmap="Reds", shade=True, bw_adjust=1)
    plt.scatter(x=medlen, y=medphred, color='black', label='({}),({})'.format(medlen,medphred))
    #plt.title('Basecalled reads length vs reads PHRED quality ')
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel("Basecalled length", size=14)
    plt.ylabel("PHRED quality scores", size=14)
    plt.xscale('log')
    fig = plt.gcf()
    p5_encoded_svg = fig_to_svg(fig)
    p5_html_svg = '<dir id="p2"><h3>Basecalled reads length vs reads PHRED quality</h3>{}<h3>Pass reads<h3></dir>'.format(p5_encoded_svg)
    plt.show()
    plt.clf()


    #5. piechart
    plt.figure()
    tools = df["barcode"].value_counts()
    # toolsitem = str(tools.index)
    # toolsval = int(tools.values)
    #print(tools,toolsitem,toolsval)
    colors = sns.color_palette('pastel')
    labels =list(df["barcode"].unique())
    plt.pie(tools, colors = colors, autopct='%.01f%%')
    #plt.legend(labels='{}'.format(df['barcode'].tolist()), loc='upper right')
    plt.legend(labels ,loc=2, bbox_to_anchor=(0.8, 1),  borderaxespad=0.)
    plt.grid()
    # #plt.pie(tools, labels = labels, colors = colors, autopct='%.0f%%')
    fig = plt.gcf()
    p7_encoded_svg = fig_to_svg(fig)
    p7_html_svg = '<dir id="p2"><h3>Number of reads per barcode</h3>{}<h3>Pass reads<h3></dir>'.format(p7_encoded_svg)
    plt.show()
    plt.clf()


    ##################################PASS READS###############################################
    ##Kernel density estimation -density plot
    #1.Kernel density estimation- Readlenght
    '''Readlenght density'''
    sns.displot(x="Readlenght", kind="kde" ,rug=True, data =passread, height=6, aspect=9/6)
    plt.xlabel("Readlenght", size=14)
    #plt.ylabel("density", size=14)
    plt.axvline(x=passread.Readlenght.mean(),ls='--',color='red', label='Mean = {}'.format(passread.Readlenght.mean()))
    plt.axvline(x=passread.Readlenght.median(),ls='--', color='black', label='Median = {}'.format(passread.Readlenght.median()))
    plt.legend(loc='upper right')
    plt.xlabel("Readlenght", size=14)
    plt.ylabel("Density", size=14)
    plt.grid()
    fig = plt.gcf()
    pr3_encoded_svg = fig_to_svg(fig)
    pr3_html_svg = '<dir id="p2"><h3>Basecalled Readlenght</h3>{}<hr></dir>'.format(pr3_encoded_svg)
    plt.show()
    plt.clf()

    #2.Kernel density estimation- MeanPhred_quality
    '''MeanPhred_quality density'''
    sns.displot(passread, x="Phred_quality", rug=True, kind="kde",height=6, aspect=9/6)
    plt.xlabel("PHRED quality scores", size=14)
    plt.ylabel("Density", size=14)
    plt.axvline(x=passread.MeanPhred_quality.mean(),ls='--',color='red', label='Mean = {}'.format(passread.MeanPhred_quality.mean()) )
    plt.axvline(x=passread.Phred_quality.median(),ls='--',color='black', label='Median = {}'.format(passread.Phred_quality.median()))
    #plt.title(' MeanPhred_quality ')
    plt.legend(loc='upper right')
    plt.grid()
    fig = plt.gcf()
    pr4_encoded_svg = fig_to_svg(fig)
    pr4_html_svg = '<dir id="p2"><h3>Basecalled reads PHRED quality</h3>{} <hr></dir>'.format(pr4_encoded_svg)
    plt.show()
    plt.clf()

    ## readlen&meanQ
    plt.figure()
    sns.kdeplot(x=passread.Readlenght, y=passread.MeanPhred_quality, cbar=True, cmap="Reds", shade=True, bw_adjust=1,)
    plt.scatter(x=passread.Readlenght.median(), y=passread.Phred_quality.median(), color='black', label='({}),({})'.format(passread.Readlenght.median(), passread.Phred_quality.median()))
    #plt.title('Basecalled reads length vs reads PHRED quality ')
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel("Basecalled length", size=14)
    plt.ylabel("PHRED quality scores", size=14)
    plt.xscale('log')
    fig = plt.gcf()
    pr5_encoded_svg = fig_to_svg(fig)
    pr5_html_svg = '<dir id="p2"><h3>Basecalled reads length vs reads PHRED quality</h3>{}<hr></dir>'.format(pr5_encoded_svg)
    plt.show()
    plt.clf()


    ###################pie
    plt.figure()
    tools2 = passread["barcode"].value_counts()
    # toolsitem = str(tools.index)
    # toolsval = int(tools.values)
    #print(tools,toolsitem,toolsval)
    colors = sns.color_palette('pastel')
    labels =list(df["barcode"].unique())
    plt.pie(tools2, colors = colors, autopct='%.01f%%')
    plt.legend(labels, bbox_to_anchor=(0.8, 1),  loc=2, borderaxespad=0.)
    plt.grid()
    fig = plt.gcf()
    pr7_encoded_svg = fig_to_svg(fig)
    pr7_html_svg = '<dir id="p2"><h3>Number of reads per barcode</h3>{}<hr></dir>'.format(pr7_encoded_svg)
    plt.show()
    plt.clf()


    ##############################HTML FILE#################################
    import webbrowser
    results = table1_html + table2_html+ p3_html_svg +pr3_html_svg+ p4_html_svg + pr4_html_svg +p5_html_svg +pr5_html_svg + p7_html_svg+ pr7_html_svg
    
    html = """
    <!doctype html>
    <html lang="en">
    <head>
        <!-- Required meta tags -->
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <meta name="description" content="Nanopore sequencing report">
        
        <!-- Bootstrap CSS -->
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-1BmE4kWBq78iYhFldvKuhfTAU6auU8tT94WrHftjDbrCEXSU1oBoqyl2QvZ6jIW3" crossorigin="anonymous">
        
    <title>Nanopore sequencing report</title>
    
        <p>This report was generated by group3 using FASTQ data from Nanopore sequencing data</p>
    
    </head>
    <body>
        <center><h1>Nanopore sequencing report</h1><center> 
    <hr>  
        
    <center>{results}<center> 
    
    <hr>  
        
    <div class="custom-select" style="width:200px;">
    <select>
        <option value="0">Select filter:</option>
        <option value="1">all reads</option>
        
        <option value="1">passed reads</option>
        </select>
    </div>
    
        
    </body>
    </html>
    """


    ########################## Write HTML String to file.html
    import webbrowser
    with open("report.html", "w" , encoding="utf-8") as file:
        file.write(html.format(results=results))
        #file.write(html.format(results=resultspass))
        file.close()
        webbrowser.open_new_tab('report.html')