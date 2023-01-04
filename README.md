Software for the interpretation of complex DNA profiles

Installation for Windows: Open R (v4.2.x) and write
install.packages('http://euroformix.com/sites/default/files/euroformix_4.0.2.zip',repos=NULL,type='win.binary')
install.packages(c('gWidgets2tcltk','cubature','XML','curl','plotly'))

Also possible to download zip file from https://github.com/oyvble/euroformix/releases (compiled for windows only)

Info about suggested R-packages (all are optional to install):\
gWidgets2tcltk: Used for R GUI\
plotly: Used to vizualiaze DNA-profiles\
cubature: Used for numerical calculation of Bayes Factor\
XML/curl: Used to read allele frequencies from Strider population databases
