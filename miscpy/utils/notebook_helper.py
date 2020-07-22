import os
from IPython.display import display, Markdown

def run_notebook_helper():
    name = 'latex_preamble.tex'
    for path in os.sys.path:
        for root, dirs, files in os.walk(path):
            if name in files:
                latex_fullpath = os.path.join(root, name)

    f = open (latex_fullpath,'r')
    txt = f.read()
    f.close()
    
    #get_ipython().set_next_input('%load $latex_fullpath')
    display(Markdown('Loading from %s\n'%(latex_fullpath)+txt))
