import papermill as pm
import sys


print("Excuting nb for res = ", 3)
pm.execute_notebook(
    'Evolution_ScalarField.ipynb',
    'Evolution_ScalarField.ipynb',
    parameters = dict(m=3,A = 0.01)
)


"""
for i in range(1, len(sys.argv)):
    print('\nResolution:', sys.argv[i],'\n\n')
"""

"""for i in range(1,3+1):
    print("Excuting nb for res = ", i)
    pm.execute_notebook(
        'Evolution_ScalarField.ipynb',
        'Evolution_ScalarField.ipynb',
        parameters = dict(m=i,A = 0.01)
    )
   
"""

"""pm.execute_notebook(
    'Evolution_ScalarField.ipynb',
    'Evolution_ScalarField.ipynb',
    parameters = dict(m=int(sys.argv[i]))
    )"""