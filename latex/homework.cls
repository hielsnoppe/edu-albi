%*******************************************************************************
%  homework.cls
%  LaTeX-Dokumentenklasse fuer Hausaufgaben
%*******************************************************************************

\NeedsTeXFormat{LaTeX2e}\relax
\ProvidesClass{homework}[2010/11/19]
\LoadClass[10pt,a4paper,landscape,twocolumn]{article}

\RequirePackage[top=2cm,right=1.5cm,bottom=1.5cm,left=1.5cm]{geometry}
\RequirePackage[utf8]{inputenc}
\RequirePackage[ngerman]{babel}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{amsthm}
\RequirePackage{listings}
\RequirePackage{graphicx}
%\RequirePackage{algo}

% meta
\newcommand{\course}[1]{\renewcommand{\course}{#1}}
\newcommand{\semester}[1]{\renewcommand{\semester}{#1}}
\newcommand{\no}[1]{\renewcommand{\no}{#1}}
\newcommand{\tutorial}[1]{\renewcommand{\tutorial}{#1}}
\newcommand{\tutor}[1]{\renewcommand{\tutor}{#1}}
\renewcommand{\date}[1]{\renewcommand{\date}{#1}}
\renewcommand{\author}[1]{\renewcommand{\author}{#1}}
\renewcommand{\maketitle}{
	\fbox{\parbox{\linewidth}{
		\course, \semester\\
		\no. \"Ubung zu \date\\
		\author\\
		Tutorium: \tutorial, \tutor
	}}
}

% convenience
\newcommand{\Nat}{\mathbb{N}}
\newcommand{\Zahl}{\mathbb{Z}}
\newcommand{\Rat}{\mathbb{Q}}
\newcommand{\Real}{\mathbb{R}}
\newcommand{\LandauO}{\mathcal{O}}
\newcommand{\Landauo}{\mathcal{o}}
\newcommand{\LandauOmega}{\Omega}
\newcommand{\Landauomega}{\omega}
\newcommand{\LandauTheta}{\Theta}
\newcommand{\floor}[1]{\left\lfloor{#1}\right\rfloor}
\newcommand{\ceil}[1]{\left\lceil{#1}\right\rceil}
\newcommand{\abs}[1]{\left|#1\right|}

%\newcounter{aufgaben}
%\newcounter{punkte}
\newcommand{\aufgabe}[2]{
	\item{\textbf{#1} \quad (#2 Punkte)\\}
%	\stepcounter{aufgaben}
%	\addtocounter{punkte}{#2}
}

% settings
\pagestyle{empty}

%\lst@definelanguage[2012]{pseudo}%
\lstdefinelanguage{pseudo}{
	morekeywords = {
		abort,and,await,begin,class,co,oc,constant,%
		datatype,do,end,ensures,false,function,if,fi,in,method,mod,%
		monitor,new,not,null,nil,od,of,or,out,override,otherwise,%
		procedure,process,rem,requires,return,spawn,skip,true,type,variable,%
		when,with,xor,yield
	},%
	sensitive = f,%
	morecomment = [l]{--},%
	morecomment = [s]{`}{`},%
	morestring = [m]",%
	morestring = [m]'%
}%[keywords,comments,strings]

\lstset{
	numberstyle=\scriptsize,
	numbersep=.5em,
	stepnumber=1,
	literate={+}{{$+$}}1
	{/}{{$/$}}1
	{*}{{$*$}}1
	{=}{{$=$}}1
	{>}{{$>$ }}1
	{<}{{$<$ }}1
	{\\\\}{{\char`\\\char`\\}}1
	{->}{{$\rightarrow$}}2
	{>=}{{$\geq$}}2
	{<-}{{$\leftarrow$}}2
	{<=}{{$\leq$}}2
	{=>}{{$\Rightarrow$ }}2
	{!=}{$\neq$}1
	{:=}{{$:=$ }}2
	{|}{{$\mid$}}1
	{\{}{{$\{$}}1
	{\}}{{$\}$}}1
	{[]}{$\square$}2,
	language = {[]pseudo},
	columns = flexible,
	basicstyle = \rmfamily\small,
	tabsize = 4
}

\parindent0mm
\sloppy
