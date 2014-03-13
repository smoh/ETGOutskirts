
path_list = strsplit(!path, ':', /extract)

plotfiles=FILE_SEARCH('/usr/home/dave/myroutines/*plot*.pro')
MK_HTML_HELP, plotfiles, 'myplot.html'

; find all *.pro files
foreach path, path_list do begin
    FILE_SEARCH(STRJOIN([path, '/*.pro']))
endforeach
