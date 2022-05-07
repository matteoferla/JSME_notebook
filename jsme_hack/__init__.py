from IPython.display import display, HTML
import ctypes

# ------------------- Deal with Colab vs Jupyter -------------------
import importlib.util

IN_COLAB = get_ipython().__class__.__name__ == 'Shell'
if not IN_COLAB:
    pass
elif importlib.util.find_spec('google'):
    from google.colab import output
else:
    raise SystemError("The IPython shell thinks it's Colab, but there's no `google.colab`!")


# ------------------- Define the class -------------------

class JSMEHack:
    """
    .. code-block:: python
      smiles = 'CCCC'
      jsme = JSMEHack('CCCC')

    ``jsme.smiles`` will change based on the JSME viewport.

    It works in Colab and Jupyter, but to work in the former, the latter is hacky.
    Namely, the id of the object is passed to JS who asks the kernel to change the ``.smiles``
    attribute of the JSMEHack instance by fetching it by its id.
    Pointers are not welcome in Python as a miscast pointer can result in a segfault.
    """

    @staticmethod
    def get_by_id(address) -> 'JSMEHack':
        return ctypes.cast(address, ctypes.py_object).value

    def set_smiles(self, smiles):
        self.smiles = smiles

    def __init__(self, smiles):
        self.smiles = smiles
        if IN_COLAB:
            output.register_callback('notebook.set_smiles', lambda smiles: self.set_smiles(smiles))
        display(HTML('''
    <script type="text/javascript" language="javascript" src="https://users.ox.ac.uk/~bioc1451/jsme/jsme.nocache.js"></script>
    <script>''' +
                     f'window.smiles = "{smiles}"; window.py_obj_id = {id(self)}' +
                     '''
                     //this function will be called after the JavaScriptApplet code has been loaded.
                         function jsmeOnLoad() {
                             const params = {smiles: smiles};
                             const jsmeApplet = new JSApplet.JSME("jsme_container", "380px", "340px", params);
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
                             Jupyter.notebook.kernel.execute(`JSMEHack.get_by_id(${py_obj_id}).smiles = '${smiles}'`);
                               }
                               else {throw "Unknown environment";}
                             });
                       }
                     </script>
                     <div id="jsme_container"></div>
                     '''))