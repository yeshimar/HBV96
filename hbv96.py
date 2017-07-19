#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HBV-96 model RoutinesClass

class HBV96(HydroModel):
	"""docstring for RoutineProcess"""

	""""""
	
	# Precipitation routine HBV96
	def _precipitation(self):
	    '''
	    ==============
	    Precipitation
	    ==============

	    Precipitaiton routine of the HBV96 model.

	    If temperature is lower than self._par['ltt'], all the precipitation is considered as
	    snow. If the temperature is higher than self._par['utt'], all the precipitation is
	    considered as rainfall. In case that the temperature is between self._par['ltt'] and
	    self._par['utt'], precipitation is a linear mix of rainfall and snowfall.

	    Parameters
	    ----------
	    intab['temp'] : float
	        Measured temperature [C]
	    self._par['ltt'] : float
	        Lower temperature treshold [C]
	    self._par['utt'] : float
	        Upper temperature treshold [C]
	    prec : float 
	        Precipitation [mm]
	    self._par['rfcf'] : float
	        Rainfall corrector factor
	    self._par['sfcf'] : float
	        Snowfall corrector factor

	    Returns
	    -------
	    _rf : float
	        Rainfall [mm]
	    _sf : float
	        Snowfall [mm]
	    '''

	    if intab['temp'] <= self._par['ltt']:
	        _rf = 0.0
	        _sf = prec*self._par['sfcf']

	    elif intab['temp'] >= self._par['utt']:
	        _rf = prec*self._par['rfcf']
	        _sf = 0.0

	    else:
	        _rf = ((intab['temp']-self._par['ltt'])/(self._par['utt']-self._par['ltt'])) * prec * self._par['rfcf']
	        _sf = (1.0-((intab['temp']-self._par['ltt'])/(self._par['utt']-self._par['ltt']))) * prec * self._par['sfcf']

	    # return _snow at line 266
		
		# Snow routine HBV96	        
		def _snow(self):
		    '''
		    ====
		    Snow
		    ====
		    
		    Snow routine of the HBV-96 model.
		    
		    The snow pack consists of two states: Water Content (wc) and Snow Pack 
		    (sp). The first corresponds to the liquid part of the water in the snow,
		    while the latter corresponds to the solid part. If the temperature is 
		    higher than the melting point, the snow pack will melt and the solid snow
		    will become liquid. In the opposite case, the liquid part of the snow will
		    refreeze, and turn into solid. The water that cannot be stored by the solid
		    part of the snow pack will drain into the soil as part of infiltration.

		    Parameters
		    ----------
		    self._par['cfmax'] : float 
		        Day degree factor
		    self._par['tfac'] : float
		        Temperature correction factor
		    intab['temp'] : float 
		        Temperature [C]
		    self._par['ttm'] : float 
		        Temperature treshold for Melting [C]
		    self._par['cfr'] : float 
		        Refreezing factor
		    self._par['cwh'] : float 
		        Capacity for water holding in snow pack
		    _rf : float 
		        Rainfall [mm]
		    _sf : float 
		        Snowfall [mm]
		    intab['wc'] : float 
		        Water content in previous state [mm]
		    intab['sp'] : float 
		        snow pack in previous state [mm]

		    Returns
		    -------
		    _in : float 
		        Infiltration [mm]
		    intab['wc'] : float 
		        Water content in posterior state [mm]
		    _intab['sp'] : float 
		        Snowpack in posterior state [mm]
		    '''

		    if intab['temp'] > self._par['ttm']:

		        if self._par['cfmax']*(intab['temp']-self._par['ttm']) < intab['sp']+_sf:
		            _melt = self._par['cfmax']*(intab['temp']-self._par['ttm'])
		        else:
		            _melt = intab['sp']+_sf

		        outab['sp'] = intab['sp'] + _sf - _melt
		        outab['wc'] = intab['wc'] + _melt + _rf

		        '''
		        # Since we begin to use hashable object to update wc value, 
		        	intermedia wc here is no longer necessary.
				
		        _wc_int = intab['wc'] + _melt + _rf 	...(1)
		        intab['wc'] = intab['wc'] + _melt + _rf 	...(2)
				
				(1) and (2) are equivalent here.
		        Modifications like this one will be applied to the whole codage.
		        '''

		    else:
		        if self._par['cfr']*self._par['cfmax']*(self._par['ttm']-intab['temp']) < intab['wc']:
		            _refr = self._par['cfr']*self._par['cfmax']*(self._par['ttm'] - intab['temp'])
		        else:
		            _refr = intab['wc'] + _rf

		        outab['sp'] = intab['sp'] + _sf + _refr
		        outab['wc'] = intab['wc'] - _refr + _rf

		    if intab['wc'] > self._par['cwh']*intab['sp']:
		        _in = outab['wc']-self._par['cwh']*intab['sp']
		        outab['wc'] = self._par['cwh']*intab['sp']
		    else:
		        _in = 0.0

		    # return _soil at line 265

		    # Soil routine HBV96
			def _soil(self):
			    '''
			    ====
			    Soil
			    ====
			    
			    Soil routine of the HBV-96 model.
			    
			    The model checks for the amount of water that can infiltrate the soil, 
			    coming from the liquid precipitation and the snow pack melting. A part of 
			    the water will be stored as soil moisture, while other will become runoff, 
			    and routed to the upper zone tank.

			    Parameters
			    ----------
			    self._par['fc'] : float 
			        Filed capacity
			    self._par['beta'] : float 
			        Shape coefficient for effective precipitation separation
			    self._par['etf'] : float 
			        Total potential evapotranspiration
			    intab['temp'] : float 
			        Temperature
			    intab['tm'] : float 
			        Average long term temperature
			    self._par['e_corr'] : float 
			        Evapotranspiration corrector factor
			    self._par['lp'] : float _soil 
			        wilting point
			    self._par['tfac'] : float 
			        Time conversion factor
			    self._par['c_flux'] : float 
			        Capilar flux in the root zone
			    _in : float 
			        actual infiltration
			    intab['ep'] : float 
			        actual evapotranspiration
			    intab['sm'] : float 
			        Previous soil moisture value
			    intab['uz'] : float 
			        Previous Upper zone value

			    Returns
			    -------
			    intab['sm'] : float 
			        New value of soil moisture
			    intab['uz'] : float 
			        New value of direct runoff into upper zone
			    '''

			    qdr = max(intab['sm'] + inf - self._par['fc'], 0)
			    _in = inf - qdr
			    _r = ((intab['sm']/self._par['fc'])**self._par['beta']) * _in
			    _ep_int = (1.0 + self._par['etf']*(intab['temp'] - intab['tm']))*self._par['e_corr']*intab['ep']
			    _ea = max(_ep_int, (intab['sm']/(self._par['lp']*self._par['fc']))*_ep_int)

			    _cf = self._par['c_flux']*((self._par['fc'] - intab['sm'])/self._par['fc'])
			    outab['sm'] = max(intab['sm'] + _in - _r + _cf - _ea, 0)
			    outab['uz'] = intab['uz'] + _r - _cf

			    # return _response at line 264

			    # Response routine HBV96
				def _response(self):
				    '''
				    ========
				    Response
				    ========
				    The response routine of the HBV-96 model.
				    
				    The response routine is in charge of transforming the current values of 
				    upper and lower zone into discharge. This routine also controls the 
				    recharge of the lower zone tank (baseflow). The transformation of units 
				    also occurs in this point.
				    
				    Parameters
				    ----------
				    self._par['tfac'] : float
				        Number of hours in the time step
				    self._par['perc'] : float
				        Percolation value [mm\hr]
				    self._par['alpha'] : float
				        Response box parameter
				    self._par['k'] : float
				        Upper zone recession coefficient
				    self._par['k1'] : float 
				        Lower zone recession coefficient
				    self._par['area'] : float
				        Catchment self._par['area'] [Km2]
				    intab['lz'] : float 
				        Previous lower zone value [mm]
				    intab['uz'] : float 
				        Previous upper zone value before percolation [mm]
				    qdr : float
				        Direct runoff [mm]
				    
				    '''    
				    outab['lz'] = intab['lz'] + np.min(self._par['perc'], outab['uz'])
				    outab['uz'] = np.max(outab['uz'] - self._par['perc'], 0.0)

				    _q_0 = self._par['k']*(outab['uz']**(1.0 + self._par['alpha']))
				    _q_1 = self._par['k1']*outab['lz']

				    outab['uz'] = max(outab['uz'] - (_q_0), 0)
				    outab['lz'] = max(outab['lz'] - (_q_1), 0)

				    outab['q_new'] = self._par['area']*(_q_0 + _q_1 + qdr)/(3.6*self._par['tfac'])

				    # return _routing at line 263

				    # Routing routine HBV96
					def _routing(self):
					    """
					    To be implemented
					    """	
					    return
					return _routing
		    	return _response
		    return _soil
	    return _snow


	def _step_run(self, intab, outab):
	    '''
	    ========
	    Step run
	    ========
	    
	    Makes the calculation of next step of discharge and states
	    
	    Parameters
	    ----------
	    p : array_like [18]
	        Parameter vector, set up as:
	        [self._par['ltt'], self._par['utt'], self._par['ttm'], self._par['cfmax'], self._par['fc'], ecorr, self._par['etf'], self._par['lp'], self._par['k'], self._par['k1'], 
	        self._par['alpha'], self._par['beta'], self._par['cwh'], self._par['cfr'], self._par['c_flux'], self._par['perc'], self._par['rfcf'], self._par['sfcf']]
	    p2 : array_like [2]
	        Problem parameter vector setup as:
	        [self._par['tfac'], self._par['area']]
	    v : array_like [4]
	        Input vector setup as:
	        [prec, intab['temp'], evap, llt]
	    St : array_like [5]
	        Previous model states setup as:
	        [sp, sm, uz, lz, wc]

	    Returns
	    -------
	    q_new : float
	        Discharge [m3/s]
	    St : array_like [5]
	        Posterior model states
	    '''

	    # Call the nested 5 routines, output will be updated input hashtables
	    _precipitation()()()()()


	def _nse(q_rec, q_sim):
	    '''
	    ===
	    NSE
	    ===
	    
	    Nash-Sutcliffe efficiency. Metric for the estimation of performance of the 
	    hydrological model
	    
	    Parameters
	    ----------
	    q_rec : array_like [n]
	        Measured discharge [m3/s]
	    q_sim : array_like [n] 
	        Simulated discharge [m3/s]
	        
	    Returns
	    -------
	    f : float
	        NSE value
	    '''
	    a = np.square(np.subtract(q_rec, q_sim))
	    b = np.square(np.subtract(q_rec, np.nanmean(q_rec)))
	    if a.any < 0.0:
	        return(np.nan)
	    f = 1.0 - (np.nansum(a)/np.nansum(b))
	    return f


	def _rmse(q_rec,q_sim):
	    '''
	    ====
	    RMSE
	    ====
	    
	    Root Mean Squared Error. Metric for the estimation of performance of the 
	    hydrological model.
	    
	    Parameters
	    ----------
	    q_rec : array_like [n]
	        Measured discharge [m3/s]
	    q_sim : array_like [n] 
	        Simulated discharge [m3/s]
	        
	    Returns
	    -------
	    f : float
	        RMSE value
	    '''
	    erro = np.square(np.subtract(q_rec,q_sim))
	    if erro.any < 0:
	        return(np.nan)
	    f = np.sqrt(np.nanmean(erro))
	    return f

	def calibrate(x_0=self._config['x_0'], obj_fun=self._config['obj_fun'],
				wu=self._config['wu'], verbose=self._config['verbose'], 
	              minimise=self._config['minimise'], fun_nam=self._config['fun_nam']):
	    '''
	    =========
	    Calibrate
	    =========

	    Run the calibration of the HBV-96. The calibration is used to estimate the
	    optimal set of parameters that minimises the difference between 
	    observations and modelled discharge.
	    
	    Parameters
	    ----------
	    
	    flow : array_like [n]
	        Measurements of discharge [m3/s]
	    intab['avg_prec'] : array_like [n]
	        Average precipitation [mm/h]
	    intab['temp'] : array_like [n]
	        Average temperature [C]
	    et : array_like [n]
	        Potential Evapotranspiration [mm/h] 
	    p2 : array_like [2]
	        Problem parameter vector setup as:
	        [self._par['tfac'], self._par['area']]
	    init_st : array_like [5], optional
	        Initial model states, [sp, sm, uz, lz, wc]. If unspecified, 
	        [0.0, 30.0, 30.0, 30.0, 0.0] mm
	    ll_temp : array_like [n], optional
	        Long term average temptearature. If unspecified, calculated from intab['temp'].
	    x_0 : array_like [18], optional
	        First guess of the parameter vector. If unspecified, a random value
	        will be sampled between the boundaries of the 
	    x_lb : array_like [18], optional
	        Lower boundary of the parameter vector. If unspecified, a random value
	        will be sampled between the boundaries of the 
	    x_ub : array_like [18], optional
	        First guess of the parameter vector. If unspecified, a random value
	        will be sampled between the boundaries of the
	    obj_fun : function, optional
	        Function that takes 2 parameters, recorded and simulated discharge. If
	        unspecified, RMSE is used.
	    wu : int, optional
	        Warming up period. This accounts for the number of steps that the model
	        is run before calculating the performance function.
	    verbose : bool, optional
	        If True, displays the result of each model evaluation when performing
	        the calibration of the hydrological model.
	    tol : float, optional
	        Determines the tolerance of the solutions in the optimisaiton process.
	    minimise : bool, optional
	        If True, the optimisation corresponds to the minimisation of the 
	        objective function. If False, the optimial of the objective function is
	        maximised.
	    fun_nam : str, optional
	        Name of the objective function used in calibration. If unspecified, is
	        'RMSE'
	    
	    Returns
	    -------
	    params : array_like [18]
	        Optimal parameter set
	    
	    performance : float
	        Optimal value of the objective function
	    '''
	    self._init_simu() # Simulation init was put here in order to save instances in _simulate_with_calibration(self)
	    
	    def _cal_fun(self, par_to_optimize):
	    	self._par.update(dict(zip(_ind[:18], par_to_optimize))) # Update the parameter dictionary
	        _q_sim, _q_rec = self._simulate_with_calibration()

	        if verbose:
	            print('{0}: {1}'.format(fun_nam, perf))
	        
	        if minimise:
	            return obj_fun(_q_rec[wu:], _q_sim[wu:])
	        else:
	            return -obj_fun(_q_rec[wu:], _q_sim[wu:])

	    # Boundaries
	    x_b = zip(self.P_LB, self.P_UB)

	    # Initial guess
	    if x_0 is None:
	        # Randomly generated
	        x_0 = np.random.uniform(self.P_LB, self.P_UB)

	    # Model optimisation
	    par_cal = opt.minimize(self._cal_fun, x_0, method='L-BFGS-B', bounds=x_b,
	                          tol=self._config['tol'])
	    self._params = par_cal.x
	    self._performance = par_cal.fun

	def _simulate_with_calibration(self):
		_q_sim = [self._data[0]['q_sim'],]
		_q_rec = [self._data[0]['q_rec'],] # _q_sim and q_rec will only live inside simulation
		
		for i in xrange(0, self._config['miles']):
			self._step_run(self._data[i], self._data[i+1]) # Consider sub-hashtable i as input and (i+1) as output table			
			_q_sim.append(self._data[i+1]['q_sim'])
			_q_rec.append(self._data[i+1]['q_rec'])
		return _q_sim, _q_rec

	def _simulate_without_calibration(self):
		self._init_simu()

		for i in xrange(0, self._config['miles']):
			self._step_run(self._data[i], self._data[i+1])

	def _init_simu(self):
		self._data[0].update(self.DEF_ST)
		self._data[0].update({'q_sim': self.DEF_q0})