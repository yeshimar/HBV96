import hbv96 as mcd

md = mcd.HBV96()

md.par['area'] = 94.36

md.par['tfac'] = 12

md.config['file_path'] = '/home/niko/Documents/UNESCO-IHE/Model/HBV96/all_data.csv'

md.config['header'] = 0

md.config['seperator'] = ','

md.config['obj_fun'] = md._rmse

md.config['init_guess'] = None

md.config['fun_name'] = 'RMSE'

md.config['wu'] = 10

md.config['verbose'] = False

md.config['minimise'] = True

md.config['tol'] = 1e7

md.dictlike_data_extractor()

md.config['miles'] = len(md.data)-1

md.calibrate()