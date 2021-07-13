from model.deepsea import DeepSEA4

class Model_Register():

    def __init__(self, model_name):
        self.model_name = model_name

    def get_model(self, n_class):

        if self.model_name == "DeepSea4":
            print("DeepSea4")
            return DeepSEA4(n_class)

        else:
            print("--------No Model Retrieved!-------")
            return
