"""
Plot the sparsity pattern of a sparse matrix.
"""
import pylab
import numpy

__docformat__ = 'restructuredtext'

def spy(A, patternonly=True, **kwargs):
    """
    To plot the sparsity pattern of a sparse matrix object. The sparse matrix
    object must implement `find()`, which must return the contents of the
    matrix in coordinate format, i.e., (val,irow,jcol).

    :parameters:

        :A:  Input matrix

        :patternonly:  If `True`, only output a black and white sparsity
                       pattern. If `False`, colorize the plot according to the
                       magnitude of the nonzero elements.

    :keywords:

        :ax: A Pylab Axes instance used to plot the sparsity pattern. If none
             is given, this function returns an Axes instance.

        :title: String to be used as title for the plot. If none is given,
                the title defaults to giving the order of the matrix and the
                number of nonzero elements.

        :comment: Comment to be appended to the default plot title in case no
                  title is supplied.

        :showtitle: Show or hide plot title.
    """
    nrow, ncol = A.shape
    sym = (A.issym != 0)
    (values,irow,jcol) = A.find()
    if patternonly:
        return FastSpy(nrow, ncol, irow, jcol, sym=sym, **kwargs)
    else:
        return FastSpy(nrow, ncol, irow, jcol, sym=sym, val=values, **kwargs)


def FastSpy(nrow, ncol, irow, jcol, **kwargs):
    """
    To plot the sparsity pattern of a sparse matrix in coordinate format.

    :parameters:

        :nrow: Number of rows of the matrix.
        :ncol: Number of columns of the matrix.
        :irow: Integer Numpy array of length nnz giving the row indices of the
               nonzero elements.
        :jcol: Integer Numpy array of length nnz giving the column indices of
               the nonzero elements.

    :keywords:

        :ax: A Pylab Axes instance used to plot the sparsity pattern. If none
             is given, this function returns an Axes instance.

        :sym: Should be set to True if the matrix is symmetric and only one
              triangle is passed in (irow,jcol).

        :title: String to be used as title for the plot. If none is given,
                the title defaults to giving the order of the matrix and the
                number of nonzero elements.

        :comment: Comment to be appended to the default plot title in case no
                  title is supplied.

        :showtitle: Show or hide plot title.

        :val: Float Numpy array of length nnz giving the values of the nonzero
              elements of the matrix. If supplied, a scatter plot is produced
              with patches of size proportional to the magnitude of the
              element. This option can slow down the plot for large values of
              the number of nonzero elements.
    """
    ax = kwargs.get('ax', None)
    if ax is None:
        fig = pylab.figure()
        ax = fig.gca()
    sym = kwargs.get('sym', False)
    t = kwargs.get('title', None)
    comment = kwargs.get('comment', '')
    showtitle = kwargs.get('showtitle', True)
    val = kwargs.get('val', None)

    nnz = len(irow)

    if val is not None:
        # Produce a scatter plot
        absval = numpy.abs(numpy.asarray(val))
        mag = 100.0/numpy.max(absval)
        s = 100.0 * absval/numpy.linalg.norm(val,ord=numpy.infty)
        ax.scatter(jcol, irow, s, s, marker='o', alpha=0.75)
        if sym:
            ax.scatter(irow, jcol, s, s, marker='o', alpha=0.75)
    else:
        ms = 1
        if nrow+ncol <= 100: ms = 5
        ax.plot(jcol, irow, 'ks', markersize=ms, linestyle='None')
        if sym:
            ax.plot(irow, jcol, 'ks', markersize=ms, linestyle='None')

    ax.xaxis.set_major_locator(pylab.MaxNLocator(nbins=9,
                                                 steps=[1, 2, 5, 10],
                                                 integer=True))
    ax.yaxis.set_major_locator(pylab.MaxNLocator(nbins=9,
                                                 steps=[1, 2, 5, 10],
                                                 integer=True))
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ticklabels = ax.get_xticklabels()
    ticklabels.extend(ax.get_yticklabels())
    if max(nrow,ncol) < 1000:
        labelfontsize = 'small'
    elif max(nrow,ncol) < 100000:
        labelfontsize = 'x-small'
    else:
        labelfontsize = 'xx-small'
    for label in ticklabels:
        label.set_fontsize(labelfontsize)
    ax.set_xlim(xmin=-1, xmax=ncol)
    ax.set_ylim(ymin=nrow, ymax=-1)
    ax.set_aspect('equal')
    if nrow == ncol:
        #ax.set_aspect('equal')
        if t is None: t = 'Order %-d with %-d nonzeros' % (nrow,nnz)
    else:
        #ratio = (1.0*nrow)/ncol
        #ax.set_aspect(ratio)
        if t is None: t = 'Size (%-d,%-d) with %-d nonzeros' % (nrow,ncol,nnz)
    if t is not None: t = comment + t
    if showtitle: ax.set_title(t, fontsize='x-small')
    return ax
