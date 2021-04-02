from QCARWorkflow.qcaw import Workflow, User, Dataset

client = User(user = "User1", password = "1234").server()

print("The list of dataset\n", client.list_collections()) # This will print all the dataset the User1 has.

ds = Dataset("test", ds_type = "OptimizationDataset", client=client).setting("load") # Loading the dataset.
print (ds.status(collapse = False))



# loading another dataset

ds_2 = Dataset("test2", ds_type = "OptimizationDataset", client=client).setting("load")
print (ds_2.status(collapse = False))
















