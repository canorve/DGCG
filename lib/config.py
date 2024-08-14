
import configparser


def read_config(confile):

	# Create a ConfigParser object
    config = configparser.ConfigParser()

	# Read the configuration file
    config.read(confile)


	# Access values from the configuration file
    # general:[General]
    img = config.get('General', 'Img')
    sexcat = config.get('General', 'SexCat')
    sigimg = config.get('General', 'SigImg',fallback="none")
    maskimg = config.get('General', 'MaskImg',fallback="none")
    psfdir = config.get('General', 'PsfDir',fallback="psfs")
    magzpt = config.getfloat('Selection', 'MagZpt',fallback=25)
    plate = config.getfloat('Selection', 'PlateScale',fallback=1)
    fileout = config.get('General', 'FileOut')


    #  GalfitSetup
    fitfunc = config.get('galfitsetup', 'fitfunc',fallback='sersic').upper()
    convbox = config.getint('GalfitSetup', 'ConvBox',fallback=100)
    fitbox = config.getint('GalfitSetup', 'FitBox',fallback=6)
    kronscale = config.getfloat('GalfitSetup', 'KronScale',fallback=1.1)
    constraints = config.getboolean('GalfitSetup', 'Constraints',fallback=True)

    try:
        imax = config.getint('GalfitSetup', 'imax', fallback=100)
        if imax < 100:
            raise ValueError("imax must be greater than 100")
    except (configparser.NoOptionError, ValueError) as e:
        print(f"Invalid imax value: {e}")
        # Handle the error, e.g., set a default value
        imax= 100


    # Selection
    try:
        galclas = config.getfloat('Selection', 'GalClas',fallback=0.6)
        if galclas > 1 or galclas < 0:
            raise ValueError("galclas must be between 0 and 1")
    except (configparser.NoOptionError, ValueError) as e:
        print(f"Invalid galclas value: {e}")
        # Handle the error, e.g., set a default value
        galclas = 0.6 


    magdiff = config.getfloat('Selection', 'MagDiff',fallback=3)
    offset = config.getint('Selection', 'Offset',fallback=10)
    magmax = config.getint('Selection', 'MagMax',fallback=17)
    maxfit = config.getint('Selection', 'MaxFit',fallback=8)
    flagsex = config.getint('Selection', 'FlagSex',fallback=4)

    try:
        nser = config.getfloat('Selection', 'NSer',fallback=1.5)
        if nser >= 10 or galclas <= 0.5:
            raise ValueError("nser must be between 0.5 and 10")
    except (configparser.NoOptionError, ValueError) as e:
        print(f"Invalid nser value: {e}")
        # Handle the error, e.g., set a default value
        nser = 1.5 



    satregionscale = config.getint('Selection', 'SatRegionScale',fallback=2)


    # Miscellaneous
    erase = config.getboolean('Miscellaneous', 'Erase',fallback=False)
    split = config.getint('Selection', 'Miscellaneous',fallback=int)
    colpar = config.get('Selection', 'ColPar',fallback="none")
    overwrite = config.getboolean('Miscellaneous', 'Overwrite',fallback=True)


    #log_level = config.get('Settings', 'log_level', fallback='INFO').upper()
    valid_fitfuncs = {'SERSIC', 'BD', 'MGE', 'ALL'}
    if fitfunc not in valid_fitfuncs:
        print(f"Invalid FitFunc: {fitfunc}. Defaulting to SERSIC.")
        fitfunc = 'SERSIC'



	# Return a dictionary with the retrieved values
    config_values = {
        #general
		'img': img,
		'sexcat': sexcat,
		'fileout': fileout,
		'fitfunc': fitfunc,
		'convbox': convbox,
		'sigimg': sigimg,
		'maskimg': maskimg,
        'psfdir': psfdir,
        'magzpt': magzpt,
        'plate': plate,

        #galfitsetup
		'fitfunc': fitfunc,
		'fitbox': fitbox,
        'convbox': convbox,
        'kronscale': kronscale,
        'constraints': constraints,
        'imax': imax,


        #selection
		'galclas': galclas, 
		'magdiff': magdiff, 
		'offset': offset, 
		'magmax': magmax, 
		'maxfit': maxfit, 
		'flagsex': flagsex, 
		'nser': nser, 
		'satregionscale': satregionscale, 


        #miscellaneous
        'erase': erase,
        'split': split,
        'colpar': colpar,
        'overwrite': overwrite,


    }


    return config_values

'''
if __name__ == "__main__":
	# Call the function to read the configuration file
    config_data = read_config()
	# Print the retrieved values
    print("image :", config_data['img']) 
    print("sextractor/clusex catalog: ", config_data['sexcat'])
    print("file output: ", config_data['fileout'])
    print("convultion box: ", config_data['convbox'])
    print("maxfit ", config_data['maxfit'])
    print("fitting box: ", config_data['fitbox'])
    print("galaxy classification: ", config_data['galclas'])
'''


