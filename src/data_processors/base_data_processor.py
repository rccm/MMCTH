class BaseDataProcessor:
    def __init__(self, data_folder, sensor):
        self.data_folder = data_folder
        self.sensor = sensor

    def read_data(self, file_name):
        # Common data reading logic
        pass

    def process_data(self):
        # Common data processing logic
        pass