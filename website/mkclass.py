import subprocess
import datetime
import os

LIBRARY = 'libnor36'


class Mkclass:

	def __init__(self, spectrum_input_file):
		self.input_spectrum_name = os.path.basename(spectrum_input_file)
		self.input_spectrum = os.path.realpath(spectrum_input_file)
		self.file_dir = os.path.dirname(os.path.realpath(spectrum_input_file))
		os.chdir(self.file_dir)
		self.prepared_spectrum = self.prepare_spectrum()
		self.is_classified = False
		print("Mk initialized")

	def prepare_spectrum(self, lambda_start="3800.0", lambda_end="5600.0", step = "0.2"):
		'''
			Returns path to the smoothed spectrum!
			Called via init, so usually no need to call this function
		'''
		#may cause problem with multiple threads?
		spectrum_input = self.input_spectrum_name

		temp_file = 'temp' + str(datetime.datetime.now())
		truncate_method = 'srebin0'
		print('Command called:', truncate_method, spectrum_input, temp_file,lambda_start, lambda_end)
		subprocess.call([truncate_method, spectrum_input, temp_file , lambda_start, lambda_end, step])
		smooth_method = 'smooth2'
		input_spacing = "0.2"
		smoothed_sp_name =  spectrum_input +'_smoothed'
		print("Command called:",smooth_method,temp_file,smoothed_sp_name,input_spacing)
		subprocess.call([smooth_method, temp_file,smoothed_sp_name, input_spacing, "3.46", "1.0"])
		print('spectrum is prepared!')
		return smoothed_sp_name


	def classify(self, library = LIBRARY, norm_type = "2", n_iter = "3"):
		'''
			Run the actual mkclass procedure on a mkclass instance
		'''
		#os.chdir(self.file_dir)
		spectrum_name = self.prepared_spectrum
		out = spectrum_name + "_out.txt"
		log = spectrum_name + "_log.txt"
		mkclass_method = "mkclass"
		print("Command called: ", mkclass_method, spectrum_name, library, out, log)
		subprocess.call([mkclass_method, spectrum_name, library , out, log, norm_type, n_iter])
		print('classified')
		self.is_classified = True


	def get_match_path(self):
		if self.is_classified == False:
			self.classify()
		return self.prepared_spectrum+ '.mat'

	def get_corrected_spec_path(self):
		'''
			Not sure how valid that is!
		'''
		if self.is_classified == False:
			self.classify()
		return self.prepared_spectrum + '.cor'

	def get_log_path(self):
		if self.is_classified == False:
			self.classify()
		return self.prepared_spectrum + '_log.txt'

	def get_out_path(self):
		if self.is_classified == False:
			self.classify()
		return self.prepared_spectrum + '_out.txt'


	def get_parsed_result(self):
		if self.is_classified == False:
			self.classify()
		with open(self.get_out_path(),'r') as f:
			raw_result = f.readline()
		raw_result = raw_result.split('|')
		result = {}
		result['class'] = raw_result[1]
		result['quality'] = raw_result[2]
		return result

def test_mk():
	test = Mkclass('files_for_mkclass/example_input.txt')
	test.classify()

if __name__ == '__main__':
	test = Mkclass('files_for_mkclass/example_input.txt')#DONT really need it 
	test.classify()
	print(test.get_out_path())
	print(test.get_parsed_result())