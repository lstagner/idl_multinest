PRO example_multinest,multinest_str

	LOG_LIKELIHOOD_FUNC = 'example_log_likelihood'
	PRIOR_FUNC = 'example_prior'
	NUM_PARAMS=2
	LIVE_NUM=300
	SAMPLE_NUM=100
	EXPAND=1.1d
	TOL=.5d
	
	multinest_str=multinest(LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,NUM_PARAMS,NUM=NUM,SAMPLE_NUM=SAMPLE_NUM,TOL=TOL,EXPAND=EXPAND,/plot)
	
END