# JSME notebook 
JSME molecular editor in a Jupyter or Colab notebook.
These two work in a very different manner and this is very inelegant, but works.

> This is not expected to be stable, but works as an interim solution.
> 
> The following discussion talks about Jupyter vs. Colab differences: 
> for more see [my blog post about it](https://blog.matteoferla.com/2022/05/js-in-colab.html)

This should be a widget (see below), but it is not for a series of reasons.

As a result I am not going to pip release it:

```bash
pip install git+https://github.com/matteoferla/JSME_notebook.git
```

[![jsme in googlecolab demo](https://img.shields.io/badge/colabs-demo.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/JSME_notebook_hack/blob/main/demo.ipynb)

I think it may be possible to make `comms` to work with it, thus removing the weirdness.
Were it done, it would no longer be that hacky.

## Example

```python
from jsme_notebook import JSMEHack

jsme = JSMEHack('CCCC')  # outputs to display...
None
```
Changes to the molecule will be available in `jsme.smiles` in subsequent cells.

For it to work in both Colab and in a Jupyter notebook,
the instance in Jupyter relies on pointers —so reassignment of the object name or deletion of the object 
will cause a segfault when the view is altered.

As an extra, there's an alternative way to use it, i.e. via Rdkit:

```python
from jsme_notebook.rdkit import JSMERDKit

mol: Chem.Mol = ...
jsme = JSMERDKit(mol)
```
## JSME in a notebook

[JSME Molecule Editor](https://jsme-editor.github.io/) is a great tool for creating and editing molecules.

There are a few Python implementations that are available for this. The other main option is Marvin for example.
The latter however requires a server to be running and return a valid SMILES string.
JSME is 100% client side.

As a result JSME is perfectly suitable for a Jupyter notebook.
Despite being an obvious candiate, there are three options for using JSME in a notebook, but only
one worked for me, but is not a package.

* [trident-chemwidgets](https://github.com/tridentbio/trident-chemwidgets) - raises a JS error
* [panel-chemistry](https://github.com/MarcSkovMadsen/panel-chemistry) - is for the Panel framework
* [JSME_ipywidget](https://github.com/lithium0003/JSME_ipywidget) - a very nice jupyter notebook
* [rdkit_ipynb_tools](https://github.com/apahl/rdkit_ipynb_tools) - I only spotted this afterwards (Jupyter notebook only)

The latter shows it is possible, but is not a module.
The first raises a JS error and has only 6 watchers, despite having documentation.
As a result I decided to make it into a package, but it turned out to be a rabbit hole of madness.

## Proper widget?

JSME_ipywidget is very interesting as it does not make a widget via the bizantine widget cookiecutter,
but circumvents it by getting `DOMWidgetView` from `'@jupyter-widgets/base'` by RequireJS
and using the fact that Jupyter serves all files within the workspace as they are listed.
Colab does not use RequireJS nor does it work for me when I load it nor does Colab serve the files 
like Jupyter nor can I get `@jupyter-widgets/base` to import.

I tried making it JSME into a proper widget, but the fact that JSME loads weirdly it did not work.
Namely:
* it looks for `jsmeOnLoad()` upon loading
* it does not export anything except pollute the namespace with `JSApplet`
* Something is weird with the npm package —and there's no working CDN

There are two packages in NPM, JSME and JSME-editor. The latter is available from 
https://cdn.skypack.dev/jsme-editor but does not work.

This explains why there is not a plethora of JSME widgets out there.

## CDN

I uploaded the JSME folder from `JSME_ipywidget` to my public html folder of my university
with a `.htaccess` file in it with `Header add Access-Control-Allow-Origin "*"`:
https://users.ox.ac.uk/~bioc1451/jsme/jsme.nocache.js 

This solves the CDN issue.

## Colab

The DOMWidgetView widget from '@jupyter-widgets/base' approach does not work for Colab
as mentioned.

As a result Colab's system is required. 

In Python

```python
from google.colab import output

def python_fun(*args):
    pass

output.register_callback('handle_for_python_fun', python_fun)
```

In Javascript

```javascript
google.colab.kernel.invokeFunction('handle_for_python_fun', [arg1, arg2])
```

Whereas in Jupyter the command for JS --> Python is in JS:

```javascript
Jupyter.notebook.kernel.execute(`python_fun(${arg1}, ${arg2})`);
```
There is the comms option, but I am not expecting much dialog,
back from Python as the JSME editor will be bound by previous cell.

The Colab handle thing means that I can use a bound method, which is fun.
With Jupyer I need to be naughty and use a pointer: i.e. deleting the faux-widget instance,
will cause a segfault.