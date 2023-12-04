from IPython.display import display, HTML
import ctypes
from warnings import warn
from typing import Optional

# ------------------- Deal with Colab vs Jupyter -------------------
import importlib.util

IN_COLAB = get_ipython().__class__.__name__ == 'Shell'
if not IN_COLAB:
    pass
elif importlib.util.find_spec('google'):
    from google.colab import output  # noqa
else:
    raise SystemError("The IPython shell thinks it's Colab, but there's no `google.colab`!")


# ------------------- Define the class -------------------

class JSMENotebook:
    """
    .. code-block:: python
      jsme = JSMENotebook('CCCC')

    ``jsme.smiles`` will change based on the JSME viewport.

    It works in Colab and Jupyter, but to work in the former, the latter is hacky.
    Namely, the id of the object is passed to JS who asks the kernel to change the ``.smiles``
    attribute of the JSMEHack instance by fetching it by its id.
    Pointers are not welcome in Python as a miscast pointer can result in a segfault.
    """
    default_cdn_url = 'https://users.ox.ac.uk/~bioc1451/jsme/jsme.nocache.js'

    @staticmethod
    def get_by_id(address) -> 'JSMEHack':
        return ctypes.cast(address, ctypes.py_object).value

    def _set_smiles(self, smiles):
        self.smiles = smiles

    def __init__(self,
                 smiles: Optional[str] = None,
                 cdn_url: Optional[str] = None,
                 ):
        if smiles is None:
            smiles = ''
        self.smiles = smiles
        if cdn_url is None:
            cdn_url = self.default_cdn_url
        if IN_COLAB:
            output.register_callback('notebook.set_smiles', lambda new_smiles: self._set_smiles(new_smiles))
        self.container_id = f'jsme_container_{hash(self)}'
        display(HTML(f'<script type="text/javascript" language="javascript" src="{cdn_url}"></script>' +
                     '<script>' +
                     f'window.smiles = "{smiles}";\n window.py_obj_id = {id(self)};\n' +
                     f'window.container_id = "{self.container_id}";\n' +
                     f'window.py_class_name = "{self.__class__.__name__}";\n' +
                     '''
                     //this function will be called after the JavaScriptApplet code has been loaded.
                         function jsmeOnLoad() {
                             const params = {smiles: smiles || undefined};
                             const jsmeApplet = new JSApplet.JSME(window.container_id, "380px", "340px", params);
                             window.jsmeApplet = jsmeApplet;
                             jsmeApplet.setCallBack("AfterStructureModified", async (jsme) => {
                               smiles = jsmeApplet.smiles();
                               // in Jupyter this is called output_area. Colab #output-area
                               // let el = document.createElement('div');
                               // el.textContent = smiles
                               // document.querySelector("#output_area,#output-area").appendChild(el);
                               if (window.google !== undefined) {
                                 await google.colab.kernel.invokeFunction('notebook.set_smiles', [smiles]);
                              } else if (window.Jupyter !== undefined) {
                              console.log(`JSMEHack.get_by_id(${py_obj_id}).smiles = '${smiles}'`);
                             Jupyter.notebook.kernel.execute(`${py_class_name}.get_by_id(${py_obj_id}).smiles = '${smiles}'`);
                               }
                               else {throw "Unknown environment";}
                             });
                       }
                     </script>'''+
                     f'<div id="{self.container_id}"></div>'
                     ))

class JSMEHack(JSMENotebook):

    def __init__(self, *args, **kwargs):
        warn('JSMEHack is now JSMENotebook', DeprecationWarning)
        super().__init__(*args, **kwargs)