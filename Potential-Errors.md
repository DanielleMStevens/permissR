## Potential Errors

The plots use Arial font, but someimes they dos not load into R properly and return this error. 
```
Error in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  : 
  invalid font type
In addition: There were 50 or more warnings (use warnings() to see the first 50)
```
This essentially means you fonts are not properly loaded. To fix this, run the following code:
```
extrafont::font_import()

or

extrafont::font_import(path = "/usr/share/fonts")
```
The propt will ask you if you want to install fonts into R, and you select 'y' and enter. Restart R and then rerun. 
