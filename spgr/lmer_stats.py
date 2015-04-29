from numpy import diag, ones, r_

from rpy2 import robjects
from rpy2.robjects.numpy2ri import activate

activate()
lme4 = robjects.packages.importr('lme4')
multcomp = robjects.packages.importr('multcomp')


def lmer(df_raw, formula='value ~ 0 + region + (1|subj)', adjust='fdr'):

    single_regions = sorted(set(df_raw['region']))

    formula = robjects.Formula(formula)
    adjustment = multcomp.adjusted('fdr')

    contr = create_contrasts(single_regions)

    lm1 = lme4.lmer(formula, data=DataFrame(df_raw))
    comps = multcomp.glht(lm1, contr)
    summary = multcomp.summary_glht(comps, test=adjustment)

    return _get_coef_pvalue(summary)


def create_contrasts(regions):

    n_regions = len(regions)

    x = diag(ones(n_regions))
    x[x == 0] = -1 / (n_regions - 1)
    x = r_[ones((1, n_regions)) * 1 / n_regions, x]

    fx = robjects.Matrix(x)
    rownames = ['intercept', ] + regions
    fx.rownames = robjects.StrVector(rownames)

    return fx


def _get_coef_pvalue(summary):
    coef_r = summary.rx2('test').rx2('coefficients')
    coefficients = dict(zip(coef_r.names, coef_r))
    coefficients = _add_intercept(coefficients)

    pvalues_r = summary.rx2('test').rx2('pvalues')
    pvalues = dict(zip(pvalues_r.names, pvalues_r))
    pvalues.pop('intercept')
    return coefficients, pvalues


def _add_intercept(coeff):
    intercept = coeff.pop('intercept')

    for k, v in coeff.items():
        coeff[k] = intercept + v
    return coeff


class DataFrame(robjects.DataFrame):
    def __init__(self, d):
        d_conv = {}
        for key, values in d.items():
            if isinstance(values[0], str):
                d_conv[key] = robjects.StrVector(values)
            elif isinstance(values[0], float) or isinstance(values[0], int):
                # very unlucky if first value is NaN or numpy
                d_conv[key] = robjects.FloatVector(values)

        super().__init__(d_conv)
