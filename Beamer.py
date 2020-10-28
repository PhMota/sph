from __future__ import print_function
import subprocess
import os

preamble = r'''
\documentclass{beamer}
%%
%% Choose how your presentation looks.
%%
%% For more themes, color themes and font themes, see:
%% http://deic.uab.es/~iblanes/beamer_gallery/index_by_theme.html
%%
\mode<presentation>
{
  \usetheme{default}      %% or try Darmstadt, Madrid, Warsaw, ...
  \usecolortheme{default} %% or try albatross, beaver, crane, ...
  \usefonttheme{default}  %% or try serif, structurebold, ...
  \setbeamertemplate{navigation symbols}{}
  \setbeamertemplate{caption}[numbered]
  \setbeamertemplate{footline}[frame number]
} 

\usepackage[english]{babel}
\usepackage[utf8x]{inputenc}
\usepackage{grffile}
\usepackage{hyperref}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{nath}
\delimgrowth=1
\setlength{\mathindent}{-1in}

\usefonttheme{professionalfonts}

\lstset{
language=Python,
breakatwhitespace=true,
breaklines=true,
basicstyle=\sffamily\tiny,
keywordstyle=\bf\color{blue},
showtabs=true,
tabsize=2,
%% digitstyle=\color{red},
%% otherstyle=\color{red},
numberstyle=\color{red},
columns=fullflexible,
%% commentstyle=\sf\tiny\color{gray},
%% stringstyle=\color{red},
showstringspaces=false,
morekeywords={as},
%% emph={(,)},
%% emphstyle=\color{blue},
numbers=left
%% identifierstyle=\color{blue}
}
\lstset{literate={-}{-}1}
\title[%s]{%s}
\author{Philipe Mota}
\institute{CBPF}
\date{\today}

\begin{document}

\begin{frame}
  \titlepage
\end{frame}
'''

B = '\n'

def pdflatex( basename ):
    print( 'pdflatex %s.tex' % basename )
    subprocess.check_output(['pdflatex', '{}.tex'.format(basename)])

#class Frame:
    

class Beamer:
    def writelines( self, s, mode = 'a' ):
        open( '%s.tex' % (self.fname+'/'+self.fname), mode ).writelines( s )
    
    def __init__(self, fname='calculations', title='simulation tools'):
        self.fname = fname
        self.basename = '{0}/{0}'.format(fname)
        self.writelines( preamble %( title, title ), mode = 'w' )
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.writelines( [r'\end{document}'] )
        pdflatex( self.basename )
    
    def par( self, s ):
        self.writelines( [s + '\n'] )

    equation_evironment = 'equation'

    def section( self, text ):
        self.writelines( [ r'\section{%s}' % text] )
    
    def environment( self, s, env = 'equation' ):
        self.writelines( [r'\begin{%s}' % env] )
        self.writelines( s )
        self.writelines( [r'\end{%s}' % env] )

    def eq( self, lhs, rhs = None, mathmode = 'equation' ):
        #print 'eq', lhs
        r = [r'\begin{%s}'%mathmode, lhs]
        if not rhs is None:
            if type(rhs) is list or type(rhs) is tuple:
                r += ['=', r'\\ '.join(rhs)]
            else:
                r += ['=', rhs]
        r += [r'\end{%s}'%mathmode]
        self.writelines( r )

    def math( self, lhs, rhs = None, label = None ):
        #print 'eq', label
        r = [r'\begin{equation}', latex(lhs)]
        if not rhs is None:
            r += ['=', latex(rhs)]
        if not label is None:
            r += [r'\label{%s}'%label]
        r += [r'\end{equation}']
        self.writelines( r )

    def matheval( self, s ):
        #print 'eq', label
        self.writelines( [r'\begin{equation}', latex(eval(m)), r'\label{%s}'%m, r'\end{equation}'] )

    def frame( self, title, *args ):
        s = r'\begin{frame}[fragile]{%s}' % title + B*3
        for arg in args:
            s += arg #+ B*2
        s += B + r'\end{frame}' + B*3
        self.writelines( s )

    def code( self, c=None, language=None, file=None, func=None ):
        if file:
            if not os.path.exists('%s/%s' % (self.fname,file)):
                print( 'running function' )
                func()
            return r'\lstinputlisting{{{}}}'.format(self.fname+'/'+file)+B
        s = r'\begin{lstlisting}[language=%s]' % (language) + B
        s += c
        s += r'\end{lstlisting}' + B
        return s
    
    def column( self, *args, **kwargs ):
        s = ''
        widths = kwargs['widths']
        s += r'\begin{columns}'+B
        for arg, width in zip(args, widths):
            s += r'\column{%s\paperwidth}'%( width if width>0 else -sum( widths ) )+B
            s += arg+B
        s += r'\end{columns}'+B
        return s
    
    def center(self, *lines ):
        s = r'\begin{center}' + B
        for line in lines:
            s += line
        s += r'\end{center}' + B
        return s
        
    def set_func( self, func ):
        self.func = func
        
    def check_file( self, fname ):
        if not os.path.exists(fname):
            print('file not found', fname)
            print( 'running function' )
            self.func()
            if not os.path.exists(fname):
                print('file not generated', fname)
                exit()
        
    def itemize( self, *items ):
        s = r'\begin{itemize}' + B
        for item in items:
            s += '\item ' + item + B
        s += r'\end{itemize}' + B
        return s
        
    
    def figure( self, path, width=None, height=None, scale=None, s='', frame=False, func=None, center=False ):
        try:
            fname = path
            self.check_file(fname)
        except:
            fname = self.fname+'/'+path
            self.check_file(fname)
        if center:
            s += r'\begin{center}' + B
        if width is not None:
            s += r'\includegraphics[width={width}\columnwidth]{{{fname}}}'.format(width=width, fname=fname) + B 
        elif height is not None:
            s += r'\includegraphics[height={height}\paperheight]{{{fname}}}'.format(height=height, fname=fname) + B
        elif scale is not None:
            s += r'\includegraphics[scale={scale}]{{{fname}}}'.format(scale=scale, fname=fname) + B
        else:
            raise Exception('missing heigh or width')
        if center:
            s += r'\end{center}' + B
        return s
        
    def table( self, file, fontsize=12, spacing=15, func=None, divide=1 ):
        s = ''
        fname = '%s/%s' % (self.fname,file)
        self.check_file(fname)
        
        for d in range(divide):
            s += r'{{\fontsize{{{}}}{{{}}}\selectfont'.format(fontsize,spacing)+B
            s += r'{\setlength{\tabcolsep}{0.2em}'+B
            for line in open(fname, 'r'):
                args = line.strip('\n').split(',')
                args = args[ d*len(args)/divide:(d+1)*len(args)/divide ]
                if '#' in line:
                    args[0] = args[0].strip('#')
                    s += r'\begin{{tabular}}{{ {} }}'.format('r'*len(args)) +B
                s += ' & '.join( args ) + r'\\' + B
            s += r'\end{tabular}'+B
            s += r'}'+B
            s += r'}'+B
        return s
    
    def tabular( self, matrix, align=None, s='' ):
        if align in ['c','r','l']:
            align = [align]*len(matrix[0])
        s += r'{\setlength{\tabcolsep}{0em}'+B
        s += r'\begin{center}' + B
        s += r'\begin{{tabular}}{{ {} }}'.format( ''.join(align) ) + B
        s += '\\\\ \n'.join( [ ' & '.join( line ) for line in matrix ] )
        s += r'\end{tabular}' + B
        s += r'\end{center}' + B
        s += r'}' + B
        #print( s )
        return s
        
def openBeamer( fname, title ): return Beamer( fname, title )

#def math( expr ):
    #return '\n'.join( ['$$', expr, '$$'] )

#def __invert__(a):
    #return 'test'

class MathGen:
    def __or__(self, a):
        return '\n'.join( ['$$', str(a), '$$'] )

#Math = MathGen()
math_counter = 0
def Math( *a ):
    global math_counter 
    math_counter += 1
    if len(a) == 1:
        a = a[0]
        if type(a) is tuple or type(a) is list:
            return '\n'.join( ['\n$$', ',\quad '.join(map(str,a)), r'.\numbered$$'+'\n'] )
        return '\n'.join( ['\n$$', str(a), r'.\numbered$$'+'\n'] )
    a = (r'\\'+'\n').join( map(str, a) )
    return '\n'.join( ['\n$$', str(a), r'.\numbered$$'+'\n'] )

def delim( a ):
    if '(' in str(a):
        if '[' in str(a):
            return Expr(r'\{{ {} \}}'.format( str(a) ))
        return Expr(r'[{}]'.format( str(a) ))
    return Expr(r'({})'.format( str(a) ))

class Expr:
    def __init__(self, s):
        self.s = s
    def __str__(self):
        return str(self.s)
    
    def __repr__(self):
        return self.s
    
    def __neg__(self):
        return Expr(r'-{}'.format(self.s))

    def __add__(self, a):
        return ExprAdd(r'{}+{}'.format(self.s, str(a)))

    def __radd__(self, a):
        return ExprAdd(r'{}+{}'.format( str(a), self.s ))

    def __sub__(self, a):
        if isinstance(a, ExprAdd):
            a = delim(a)
        return ExprAdd(r'{}-{}'.format(self.s, str(a)))

    def __rsub__(self, a):
        return ExprAdd(r'{}-{}'.format( str(a), self.s ))
    
    def __div__( self, a ):
        return frac(self.s, a)

    def __rdiv__( self, a ):
        return frac(a, self.s )
    
    def __mul__( self, a ):
        if isinstance( a, ExprAdd ):
            return a.__rmul__(self.s)
        return ExprMul(r'{}{}'.format(self.s, str(a) ))

    def __rmul__( self, a ):
        return ExprMul(r'{}{}'.format( str(a), self.s ))
    
    def __or__( self, a ):
        return Expr(r'{}|{}'.format(self.s, str(a) ))
    
    def __call__( self, *a ):
        p = ', '.join( map(str,a) )
        p = delim(p)
        return Expr(r'{}{}'.format( self.s, p ) )
        
    def __eq__( self, a ):
        #print( self.s, 'Expr.__eq__', str(a) )
        if isinstance(a, ExprEq):
            #print( self.s, 'Expr.__eq__(ExprEq)', str(a) )
            #print( ExprEq(r'{} = {}'.format( self.s, str(a) ) ) )
            return ExprEq(r'{} = {}'.format( self.s, str(a) ) )
        if type(a) is tuple or type(a) is list:
            return ExprEq(r'{} \wall = {} \return'.format( self.s, r' \\ ='.join(map(str,a)) ) )
        #if isinstance( a, Expr ):
            #print( self.s, 'Expr.__eq__(Expr)', str(a) )
            #print( ExprEq(r'{} = {}'.format( self.s, str(a) ) ) )
        #self.s = ExprEq(r'{} = {}'.format( self.s, str(a) ) )
        #print( self.s, '__eq__', a )
        #a.s = ExprEq(r'{} = {}'.format( self.s, str(a) ) )
        return ExprEq(r'{} = {}'.format( self.s, str(a) ) )

    def __bool__(self):
        print( '__bool__' )
        return True
    
    def __and__( self, a ):
        print( self.s, '__and__', str(a) )
        return Expr(self.s)

    def and_( self, a ):
        print( self.s, 'and_', str(a) )
        return Expr(self.s)
        
    def __getitem__( self, a ):
        if type(a) is tuple:
            return Expr(r'{}_{{{}}}^{{{}}}'.format( self.s, str(a[0]), str(a[1]) ) )
        return Expr(r'{}_{{ {} }}'.format( self.s, str(a) ) )
    
    def __pow__( self, a ):
        return Expr(r'{}^{{ {} }}'.format( self.s, str(a) ) )

    def __rpow__( self, a ):
        return Expr(r'{}^{{ {} }}'.format( str(a), self.s ) )

    def __sqrt__(self):
        return Expr(r'\sqrt{{ {} }}'.format( self.s ) )
    
    def __nonzero__(self):
        print( 'Expr.__nonzero__' )
        return 1
    
class ExprGen:
    def __or__( self, a ):
        return Expr(a)

class ExprEq(Expr):
    def __init__(self, s):
        #print( 'ExprEq init' )
        self.s = s

    def and_(self, a):
        print( self.s, 'ExprEq.and_', str(a) )

    def __and__(self, a):
        print( self.s, 'ExprEq.__and_', str(a) )
        return ExprEq( r'{{}}, \quad {{}}'.format(self.s, str(a) ) )

    def __eq__(self, a):
        global math_counter
        #print( self.s, 'ExprEq.__eq__', str(a), math_counter )
        #print( ExprEq( r'{} = {}'.format(self.s, str(a) ) ) )
        self.s = ExprEq( r'{} = {}'.format(self.s, str(a) ) )
        return self.s
    
    def __nonzero__(self):
        #print( '__nonzero__(ExprEq)', self.s )
        return 1
    
class ExprAdd(Expr):
    def __init__(self, s):
        self. s = s

    def __pow__( self, a ):
        return Expr(r'{}^{}'.format( delim(self.s), str(a) ) )
    
    def __mul__( self, a ):
        return ExprMul(r'{}{}'.format( delim(self.s), str(a) ) )

    def __rmul__( self, a ):
        return ExprMul(r'{}{}'.format( str(a), delim(self.s) ) )
    
    def __neg__(self):
        return Expr(r'-{}'.format(delim(self.s)))
    
    def __rsub__(self, a):
        return Expr(r'{}-{}'.format(str(a),delim(self.s)) )
    
    
class ExprMul(Expr):
    def __init__(self, s):
        self. s = s
    def __pow__( self, a ):
        return Expr(r'{}^{{ {} }}'.format( delim(self.s), str(a) ) )

def matrix( a ):
    return Expr( 
        '\n'.join([ 
            r'[\begin{matrix}',
            (r'\\' + '\n').join(
                map(
                    lambda line: ' & '.join( map(
                        str, line
                        )),
                    a)
                ),
            r'\end{matrix}]',
            ]) 
        )
        
    
mu = Expr(r'\mu ')
sigma = Expr(r'\sigma ')
alpha = Expr(r'\alpha ')
pi = Expr(r'\pi ')
xi = Expr(r'\xi ')

Gamma = Expr(r'\Gamma ')
Prod = Expr(r'\prod ')
Sum = Expr(r'\sum ')
exp = Expr(r'{\rm exp}')
e = Expr(r'{\rm e}')
log = Expr(r'{\rm log}')
erf = Expr(r'{\rm erf}')
partial = Expr(r'\partial ')

Int = Expr(r'\int')
inf = Expr(r'\infty')
oo = Expr(r'\infty')

E = ExprGen()

def cross( a, b ):
    if isinstance(a, ExprAdd):
        a = delim(a)
    if isinstance(b, ExprAdd):
        b = delim(b)
    return Expr( r'{}\times {}'.format( str(a), str(b)) )

def cdot( a, b ):
    if isinstance(a, ExprAdd):
        a = delim(a)
    if isinstance(b, ExprAdd):
        b = delim(b)
    return Expr( r'{}\cdot {}'.format( str(a), str(b)) )

def frac(a, b):
    return Expr( r'\frac{{ {} }}{{ {} }}'.format( str(a), str(b)) )

def rm(a):
    return Expr( r'{{\rm {} }}'.format(str(a)) )

def bf(a):
    return Expr( r'{{\mathbf {} }}'.format(str(a)) )

def ave(a):
    return Expr( r'\langle {} \rangle '.format(str(a)) )

def bar(a):
    return Expr( r'\bar{{ {} }}'.format(str(a)) )

def dot(a):
    return Expr( r'\dot{{ {} }}'.format(str(a)) )

def hat(a):
    return Expr( r'\hat{{ {} }}'.format(str(a)) )

def cal(a):
    return Expr( r'{{\mathcal {} }}'.format(str(a)) )

def sqrt(a):
    return Expr(r'\sqrt{{ {} }}'.format( str(a) ) )

def d(a):
    return Expr(r'{{\rm d}} {}\, '.format( str(a) ) )

def factorial(a):
    if isinstance(a, ExprAdd) or isinstance(a, ExprMul):
        Expr(r'({})!'.format( str(a) ) )
    return Expr(r'{}!'.format( str(a) ) )
