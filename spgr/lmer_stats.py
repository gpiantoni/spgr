from numpy import diag, ones, r_

from rpy2 import robjects
from rpy2.robjects.numpy2ri import activate
from rpy2.robjects.packages import importr


activate()
lme4 = importr('lme4')
multcomp = importr('multcomp')


def add_to_dataframe(df, subj, values, chan):
    """Add values for each electrode to the main frame.
    """
    for i, one_chan in enumerate(chan.chan):
        one_value = values[i]
        region = one_chan.attr['region']
        if region[:3] == 'ctx':
            df['subj'].append(subj)
            df['region'].append(region[7:])
            df['elec'].append(one_chan.label)
            df['value'].append(one_value)


def lmer(df_raw, lg, formula='value ~ 0 + region + (1|subj)', adjust='fdr',
         pvalue=0.05):
    """Compute linear mixed-effects models, using R

    Parameters
    ----------
    df_raw : dict
        dict where each key is one column
    lg : instance of logging.Logger
        logging template
    formula : str
        formula to test
    adjust : str
        adjustment ('fdr')
    pvalue : float
        threshold for p-value to report the result.

    Returns
    -------
    dict
        dictionary with coefficients
    dict
        dictionary with pvalues
    """
    single_regions = sorted(set(df_raw['region']))

    formula = robjects.Formula(formula)
    adjustment = multcomp.adjusted(adjust)

    contr = _create_contrasts(single_regions)

    lm1 = lme4.lmer(formula, data=DataFrame(df_raw))
    comps = multcomp.glht(lm1, contr)
    summary = multcomp.summary_glht(comps, test=adjustment)

    coef, intercept, pvalues = _get_coef_pvalue(summary)
    _report_values(lg, coef, pvalues, intercept, pvalue)

    return coef, pvalues


def _report_values(lg, coef, pvalues, intercept, p_threshold):
    """Report values of the LMER statistics, including the intercept value

    Parameters
    ----------
    lg : instance of logging.Logger
        logging template
    coef : dict
        dictionary with coefficients
    pvalues : dict
        dictionary with pvalues
    intercept : float
        value for the intercept
    p_threshold : float
        threshold for p-value to show values
    """
    lg.info('LMER summary')
    has_intercept = False
    for region, _ in sorted(coef.items(), key=lambda x: x[1]):
        if coef[region] > intercept and not has_intercept:
            lg.info('{:30} coef={:.3f}'.format('', intercept))
            has_intercept = True

        if pvalues[region] <= p_threshold:
            lg.info('{:30} coef={:.3f},  p-value = {:.3f}'
                    ''.format(region, coef[region], pvalues[region]))

    if not has_intercept:  # in theory, this should never happen
        lg.info('{:16} coef={:.3f}'.format('', intercept))


def _create_contrasts(regions):
    """Create contrasts

    Parameters
    ----------
    regions : list of str
        regions to create contrast matrix

    Returns
    -------
    robjects.Matrix
        matrix with the contrasts, with size (n_regions+1) X n_regions. We can
        use this to compute the intercept as well as the other contrasts (maybe
        one too many).
    """
    n_regions = len(regions)

    x = diag(ones(n_regions))
    x[x == 0] = -1 / (n_regions - 1)
    x = r_[ones((1, n_regions)) * 1 / n_regions, x]

    fx = robjects.Matrix(x)
    rownames = ['intercept', ] + regions
    fx.rownames = robjects.StrVector(rownames)

    return fx


def _get_coef_pvalue(summary):
    """Get coefficients and p-values

    Parameters
    ----------
    summary : robjects.S4
        summary from summary_glht

    Returns
    -------
    dict
        dictionary with coefficients
    float
        value for the intercept
    dict
        dictionary with pvalues
    """
    coef_r = summary.rx2('test').rx2('coefficients')
    coefficients = dict(zip(coef_r.names, coef_r))
    intercept = coefficients['intercept']
    coefficients = _add_intercept(coefficients)

    pvalues_r = summary.rx2('test').rx2('pvalues')
    pvalues = dict(zip(pvalues_r.names, pvalues_r))
    pvalues.pop('intercept')

    return coefficients, intercept, pvalues


def _add_intercept(coeff):
    """Add the intercept value back into the estimates, so that it's easier to
    interpret.

    Parameters
    ----------
    coeff : dict
        dictionary with coefficients

    Returns
    -------
    dict
        dictionary with coefficients, where the values also include the
        intercept (and intercept has been removed from dict)
    """
    intercept = coeff.pop('intercept')

    for k, v in coeff.items():
        coeff[k] = intercept + v
    return coeff


class DataFrame(robjects.DataFrame):
    """Dataframe to convert dictionary into robjections.Dataframe

    Parameters
    ----------
    d : dict
        dictionary with values as float or str
    """
    def __init__(self, d):
        d_conv = {}
        for key, values in d.items():
            if isinstance(values[0], str):
                d_conv[key] = robjects.StrVector(values)
            elif isinstance(values[0], float) or isinstance(values[0], int):
                # very unlucky if first value is NaN or numpy
                d_conv[key] = robjects.FloatVector(values)

        super().__init__(d_conv)
