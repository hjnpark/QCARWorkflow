from QCARWorkflow.qcaw import Workflow, User, Dataset

client = User(user = "User1", password = "1234").server()

print("The list of dataset\n", client.list_collections()) # This will print all the dataset the User1 has.

ds = Dataset("test", ds_type = "OptimizationDataset", client=client).setting("make") # Creating an OptimizationData set named 'test'.
wf = Workflow(ds = ds, client = client, spec_name = "spec_test", initial = "HCN.xyz") # Workflow class
wf.optimization(method = 'b3lyp', basis ='def2-svpd', subsample = 1, compute = True) # optimization function
print (ds.status(collapse = False))



# We can submit another set of molecules with a new data set

ds_2 = Dataset("test2", ds_type = "OptimizationDataset", client=client).setting("make")
Workflow(ds = ds_2, client = client, spec_name = "spec_test2", initial = "H2O_proton.xyz").optimization(method = 'b3lyp', basis = 'def2-svpd', subsample = 20, compute = True)
print (ds_2.status(collapse = False))
















