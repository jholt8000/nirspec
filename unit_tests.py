
import nirspec_wavelength_utils
import unittest
import nirspec_test_data

#TODO need to define dx, dy inputs for these tests to work

### nirspec_wavelength_utils tests ###
def test_read_oh():
    # do not need dx and, just that low_d is False
    t1 = nirspec_wavelength_utils.LineId(False, [], [])
    t1.read_OH()
    if t1.ohx = nirspec_test_data.hardcoded.ohx:
        print 'nirspec_wavelength_utils.read_OH ohx passed'
    else:
        print 'nirspec_wavelength_utils.read_OH ohx failed'

    if t1.ohy = nirspec_test_data.hardcoded.ohy:
        print 'nirspec_wavelength_utils.read_OH ohy passed'
    else:
        print 'nirspec_wavelength_utils.read_OH ohy failed'

def test_gauss_sky():
    # does not need dx, dy, or low_d, perhaps this method belongs elsewhere?
    t1 = nirspec_wavelength_utils.LineId(False, [], [])
    ohx=[100,200,220.1, 300]
    ohy=[0.1, 10, 50, 74.3]
    t1.gauss_sky(ohx, ohy, sig=0.2)
    return t1.fake_sky

def test_find_xcorr_shift(ohgs):
    t1 = nirspec_wavlelength_utils.LineId(False, dx, dy)


def test_identify(self, ohx, ohy):
    pass

def test_sanity_check(orig_pix_x, order_number_array, matched_sky_line):
    pass


### nirspec_util.nirspecBookkeeping

def test_close_logger(self):
       pass

def test_setup_logger(self):
        """
        creates logging utility, returns logger object
        """
        log = logging.getLogger(__name__)
        x = list(log.handlers)
        for i in x:
            log.removeHandler(i)
            i.flush()
            i.close()

        if '/' in self.fits_name:
            fits_nameforlog = self.fits_name.split('/')[-1]
        else:
            fits_nameforlog = self.fits_name
        if os.path.isfile(self.outpath + '/drp_' + fits_nameforlog.rstrip('.fits') + '.log'):
            if os.path.isfile(self.outpath + '/drp_' + fits_nameforlog.rstrip('.fits') + '_old.log'):
                os.remove(self.outpath + '/drp_' + fits_nameforlog.rstrip('.fits') + '_old.log')
            os.rename(self.outpath + '/drp_' + fits_nameforlog.rstrip('.fits') + '.log',
                      self.outpath + 'drp_' + fits_nameforlog.rstrip('.fits') + '_old.log')

        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

        if not logger.handlers:
            fh = logging.FileHandler(filename=self.outpath + 'drp_' + fits_nameforlog.rstrip('.fits') + '.log')
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(formatter)
            logger.addHandler(fh)

            if self.verbose:
                ch = logging.StreamHandler()
                ch.setLevel(logging.DEBUG)
                ch.setFormatter(formatter)
                logger.addHandler(ch)
            else:
                logger.propagate = False

        return logger

### nirspec_util no class
def test_make_nirspec_final_fits_and_plots(allreduceobj, order_num, sciorder, lineobj, traceobj, flatobj, sciobj):


### nirspec_util.NirspecHeader need to init a header

def test_get_data_dict(self):
    pass

def test_get_theory_order_pos(self, order_num, a_to_mu=True):
    pass

def test_sanity_check(self, flatheader):
    pass

def test_get_actual_order_num_pos(edges, theory, threshold):
    pass


#### twod_lambda_fit.py
def test_twod_lambda_fit.twodfit(dataX, dataY, dataZ, logger, lower_len_points=10., sigma_max=0.5):
    pass

### Test repo

### reduce_order

def test_reduce_order(reduction, sciobj, flatobj):
    pass


### trace_order
def __init__(self, reductionobj, flatobj, sciobj):
    self.reductionobj = reductionobj
    self.flatobj = flatobj
    self.sciobj = sciobj
    self.trace_success = False

def trace_order(self):
    pass

def test_fudge_padding(self):
    pass

def test_smooth_spectroid(self):
    pass


def test_shift_order(self):
    pass


def test_shift_order_back(self):
    pass

def test_determine_lhs_edge_pos(self, edges, theory_lhs, data_dict):
    pass