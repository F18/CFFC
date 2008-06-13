
/*----------------------------------------------------------------------------

   Description  :   C interface to gnuplot
  
    gnuplot is a freely available, command-driven graphical display tool for
    Unix. It compiles and works quite well on a number of Unix flavours as
    well as other operating systems. The following module enables sending
    display requests to gnuplot through simple C calls.
  
 ---------------------------------------------------------------------------*/

#ifndef _GNUPLOT_INCLUDED
#define _GNUPLOT_INCLUDED

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>



#define GP_MAX_TMP_FILES    64
#define GP_TMP_NAME_SIZE    512
#define GP_CMD_SIZE     	1024

/*---------------------------------------------------------------------------
                                New Types
 ---------------------------------------------------------------------------*/

#ifndef _ECLIPSE_TYPES_H_
/*
 * dpoint is convenient to store signals which have definition both on x and
 * y axis.
 */
typedef struct _DPOINT_ {
    double  x ;
    double  y ;
} dpoint ;
#endif


/*
 * This structure holds all necessary information to talk to a gnuplot
 * session.
 */

class Gnuplot_Control {
public:
    FILE *gnucmd;
    
    int nplots;
    char pstyle[32];
    
    char to_delete[GP_MAX_TMP_FILES][GP_TMP_NAME_SIZE];
    int ntmp;
    
    
    int check_X_display(int activate);
    char * gnuplot_get_program_path(char * pname);
    int  gnuplot_init(void);
    void gnuplot_close(void);
    void gnuplot_cmd(char *  cmd, ...);
    void gnuplot_setstyle(char * plot_style);
    void gnuplot_resetplot(void);
    void gnuplot_plot1d_var1(double *d,
                             int n_point,
                             char *title) ;
    void gnuplot_plot1d_var2(dpoint *d,
                             int n_points,
                             char *title) ;
    void gnuplot_plot1d_var2(double *x,
                             double *y,
                             int n_points,
                             char *title);
    void gnuplot_plot_slope(double a,
                            double b,
                            char *title) ;
    void gnuplot_plot_equation(char *equation,
                               char *title) ;
    void gnuplot_set_xlabel(char *label);
    void gnuplot_set_ylabel(char *label);
    void gnuplot_set_title(char *title);
    void gnuplot_set_title(string title);


};




/*---------------------------------------------------------------------------
 Defines
 ---------------------------------------------------------------------------*/

/* Maximal size of a gnuplot command */
#define GP_CMD_SIZE     1024
/* Maximal size of a plot title */
#define GP_TITLE_SIZE   80
/* Maximal size for an equation */
#define GP_EQ_SIZE      512


/*---------------------------------------------------------------------------
 Function codes
 ---------------------------------------------------------------------------*/




inline int Gnuplot_Control::check_X_display ( int activate )
/*-------------------------------------------------------------------------*/
/**
 @name		check_X_display
 @memo		Checks out if the DISPLAY environment variable is set.
 @param	activate int flag
 @return	int 1 if the variable is set, 0 otherwise.
 @doc
 
 This function checks out the DISPLAY environment variable to see if
 it exists. It does not check if the display is actually correctly
 configured. If you do not want to activate this check (e.g. on
 systems that do not support this kind of display mechanism), pass a
 0 integer as the activate flag. Any other value will activate it.
 */
/*--------------------------------------------------------------------------*/
{
    char *display;
    
    if ( !activate ) 
    {
        return 1;
    }
    
    display = getenv ( "DISPLAY" );
    
    if ( display == NULL )
    {
        fprintf ( stderr, "cannot find DISPLAY variable: is it set?\n");
        return 1 ;
    } 
    else
    {
        return 0 ;
    }
}

/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_get_program_path
 @memo		Find out where a command lives in your PATH.
 @param	pname Name of the program to look for.
 @return	pointer to statically allocated character string.
 @doc
 
 This is the C equivalent to the 'which' command in Unix. It parses
 out your PATH environment variable to find out where a command
 lives. The returned character string is statically allocated within
 this function, i.e. there is no need to free it. Beware that the
 contents of this string will change from one call to the next,
 though (as all static variables in a function).
 
 The input character string must be the name of a command without
 prefixing path of any kind, i.e. only the command name. The returned
 string is the path in which a command matching the same name was
 found.
 
 Examples (assuming there is a prog named 'hello' in the cwd):
 
 \begin{itemize}
 \item gnuplot_get_program_path("hello") returns "."
 \item gnuplot_get_program_path("ls") returns "/bin"
 \item gnuplot_get_program_path("csh") returns "/usr/bin"
 \item gnuplot_get_program_path("/bin/ls") returns NULL
 \end{itemize}
 
 */
/*-------------------------------------------------------------------------*/

#define MAXNAMESZ       4096
inline char * Gnuplot_Control::gnuplot_get_program_path(char * pname)
{
    int         i, j, lg;
    char    *   path;
    static char buf[MAXNAMESZ];
    
    /* Trivial case: try in CWD */
    sprintf(buf, "./%s", pname) ;
    if (access(buf, X_OK)==0) {
        sprintf(buf, ".");
        return buf ;
    }
    /* Try out in all paths given in the PATH variable */
    buf[0] = 0;
    path = getenv("PATH") ;
    if (path!=NULL) {
        for (i=0; path[i]; ) {
            for (j=i ; (path[j]) && (path[j]!=':') ; j++);
            lg = j - i;
            strncpy(buf, path + i, lg);
            if (lg == 0) buf[lg++] = '.';
            buf[lg++] = '/';
            strcpy(buf + lg, pname);
            if (access(buf, X_OK) == 0) {
                /* Found it! */
                break ;
            }
            buf[0] = 0;
            i = j;
            if (path[i] == ':') i++ ;
        }
    }
    /* If the buffer is still empty, the command was not found */
    if (buf[0] == 0) return NULL ;
    /* Otherwise truncate the command name to yield path only */
    lg = strlen(buf) - 1 ;
    while (buf[lg]!='/') {
        buf[lg]=0 ;
        lg -- ;
    }
    buf[lg] = 0;
    
    return buf ;
}
#undef MAXNAMESZ



/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_init
 @memo		Opens up a gnuplot session, ready to receive commands.
 @return	Newly allocated gnuplot control structure.
 @doc
 
 This opens up a new gnuplot session, ready for input. The struct
 controlling a gnuplot session should remain opaque and only be
 accessed through the provided functions.
 */
/*--------------------------------------------------------------------------*/

inline int Gnuplot_Control::gnuplot_init ( void )
{
    if (check_X_display(1)) return NULL ;
    
//	if (gnuplot_get_program_path("gnuplot")==NULL) {
//        fprintf(stderr, "cannot find gnuplot in your PATH");
//        return 1;
//	}
    
    /* 
     * Structure initialization:
     */
    nplots = 0 ;
    gnuplot_setstyle("points") ;
    ntmp = 0 ;
    
    gnucmd = popen("gnuplot", "w") ;
    //gnucmd = popen("/usr/local/bin/gnuplot", "w") ;
    
    if (gnucmd == NULL) {
        fprintf(stderr, "error starting gnuplot\n") ;
        return 1;
    }
    return 0;
}


/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_close
 @memo		Closes a gnuplot session previously opened by gnuplot_init()
 @param	handle Gnuplot session control handle.
 @return	void
 @doc
 
 Kills the child PID and deletes all opened temporary files.
 It is mandatory to call this function to close the handle, otherwise
 temporary files are not cleaned and child process might survive.
 
 */
/*--------------------------------------------------------------------------*/

inline void Gnuplot_Control::gnuplot_close( void )
{
    int     i ;
    if (check_X_display(1)) return ;
    if (ntmp) {
        for (i=0 ; i<ntmp ; i++) {
            remove(to_delete[i]) ;
        }
    }
    if (pclose(gnucmd) == -1) {
        fprintf(stderr, "problem closing communication to gnuplot\n") ;
        return ;
    }
    return ;
}


/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_cmd
 @memo		Sends a command to an active gnuplot session.
 @param	handle Gnuplot session control handle
 @param	cmd    Command to send, same as a printf statement.
 @return	void
 @doc
 
 This sends a string to an active gnuplot session, to be executed.
 There is strictly no way to know if the command has been
 successfully executed or not.
 The command syntax is the same as printf.
 
 Examples:
 
 \begin{itemize}
 \item gnuplot_cmd(g, "plot %d*x", 23.0);
 \item gnuplot_cmd(g, "plot %g * cos(%g * x)", 32.0, -3.0);
 \end{itemize}
 
 */
/*--------------------------------------------------------------------------*/

inline void Gnuplot_Control::gnuplot_cmd (char *  cmd, ...)
{
    va_list ap ;
    char    local_cmd[GP_CMD_SIZE];
    
    va_start(ap, cmd);
    vsprintf(local_cmd, cmd, ap);
    va_end(ap);
    
    strcat(local_cmd, "\n");
    
    fputs(local_cmd, gnucmd) ;
    fflush(gnucmd) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_setstyle
 @memo		Change the plotting style of a gnuplot session.
 @param	h Gnuplot session control handle
 @param	plot_style Plotting-style to use (character string)
 @return	void
 @doc
 
 The provided plotting style is a character string. It must be one of
 the following:
 
 \begin{itemize}
 \item {\it lines}
 \item {\it points}
 \item {\it linespoints}
 \item {\it impulses}
 \item {\it dots}
 \item {\it steps}
 \item {\it errorbars}
 \item {\it boxes}
 \item {\it boxeserrorbars}
 \end{itemize}
 */
/*--------------------------------------------------------------------------*/

inline void Gnuplot_Control::gnuplot_setstyle(char * plot_style) 
{
    if (strcmp(plot_style, "lines") &&
        strcmp(plot_style, "points") &&
        strcmp(plot_style, "linespoints") &&
        strcmp(plot_style, "impulses") &&
        strcmp(plot_style, "dots") &&
        strcmp(plot_style, "steps") &&
        strcmp(plot_style, "errorbars") &&
        strcmp(plot_style, "boxes") &&
        strcmp(plot_style, "boxerrorbars")) {
        fprintf(stderr, "warning: unknown requested style: using points\n") ;
        (void)strcpy(pstyle, "points") ;
    } else {
        (void)strcpy(pstyle, plot_style) ;
    }
    return ;
}


/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_set_xlabel
 @memo		Sets the x label of a gnuplot session.
 @param	h Gnuplot session control handle.
 @param	label Character string to use for X label.
 @return	void
 @doc
 
 Sets the x label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

inline void Gnuplot_Control::gnuplot_set_xlabel(char * label)
{
    char    cmd[GP_CMD_SIZE] ;
    
    (void)sprintf(cmd, "set xlabel \"%s\"", label) ;
    gnuplot_cmd(cmd) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_set_ylabel
 @memo		Sets the y label of a gnuplot session.
 @param	h Gnuplot session control handle.
 @param	label Character string to use for Y label.
 @return	void
 @doc
 
 Sets the y label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

inline void Gnuplot_Control::gnuplot_set_ylabel(char * label)
{
    char    cmd[GP_CMD_SIZE] ;
    
    (void)sprintf(cmd, "set ylabel \"%s\"", label) ;
    gnuplot_cmd(cmd) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_set_title
 @memo		Sets the y label of a gnuplot session.
 @param	h Gnuplot session control handle.
 @param	title Character string to use for the title.
 @return	void
 @doc
 
 Sets the y label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

inline void Gnuplot_Control::gnuplot_set_title(char * title)
{
    char    cmd[GP_CMD_SIZE] ;
    
    (void)sprintf(cmd, "set title \"%s\"", title) ;
    gnuplot_cmd(cmd) ;
    return ;
}
inline void Gnuplot_Control::gnuplot_set_title(string title)
{
    char    cmd[GP_CMD_SIZE] ;
    
    (void)sprintf(cmd, "set title \"%s\"", title.c_str()) ;
    gnuplot_cmd(cmd) ;
    return ;
}

/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_resetplot
 @memo		Resets a gnuplot session (next plot will erase previous ones).
 @param	h Gnuplot session control handle.
 @return	void
 @doc
 
 Resets a gnuplot session, i.e. the next plot will erase all previous
 ones.
 */
/*--------------------------------------------------------------------------*/

inline void Gnuplot_Control::gnuplot_resetplot(void)
{
    int     i ;
    if (ntmp) {
        for (i=0 ; i<ntmp ; i++) {
            remove(to_delete[i]) ;
        }
    }
    ntmp = 0 ;
    nplots = 0 ;
    return ;
}



/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_plot1d_var1
 @memo		Plots a 2d graph from a list of doubles.
 @param	d		Pointer to a list of doubles.
 @param	n_point	Number of doubles in the list.
 @param	title	Title of the plot.
 @return	void
 @doc
 
 Plots out a 2d graph from a list of doubles. The x-coordinate is the
 index of the double in the list, the y coordinate is the double in
 the list.
 
 Example:
 
 \begin{verbatim}
 gnuplot_ctrl    *h ;
 double          d[50] ;
 int             i ;
 
 h = gnuplot_init() ;
 for (i=0 ; i<50 ; i++) {
 d[i] = (double)(i*i) ;
 }
 gnuplot_plot1d_var1(h, d, 50, "parabola") ;
 sleep(2) ;
 gnuplot_close(h) ;
 \end{verbatim}
 
 */
/*--------------------------------------------------------------------------*/

inline void Gnuplot_Control::gnuplot_plot1d_var1(
                                          double          *   d,
                                          int                 n_point,
                                          char            *   title
)
{
    int         i ;
    FILE    *   tmp ;
    char    *   name ;
    char        cmd[GP_CMD_SIZE] ;
    char        line[GP_CMD_SIZE] ;
    
    /* can we open one more temporary file? */
    if (ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return ;
    }
    
    /* Open temporary file for output   */
    if ((name = tmpnam(NULL)) == (char*)NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }
    if ((tmp = fopen(name, "w")) == NULL) {
        fprintf(stderr, "cannot create temporary file: exiting plot") ;
        return ;
    }
    
    /* Store file name in array for future deletion */
    (void)strcpy(to_delete[ntmp], name) ;
    ntmp ++ ;
    
    /* Write data to this file  */
    for (i=0 ; i<n_point ; i++) {
        (void)fprintf(tmp, "%g\n", d[i]) ;
    }
    (void)fflush(tmp) ;
    (void)fclose(tmp) ;
    
    /* Command to be sent to gnuplot    */
    if (nplots > 0) {
        (void)strcpy(cmd, "replot") ;
    } else {
        (void)strcpy(cmd, "plot") ;
    }
    
    if (title == NULL) {
        (void)sprintf(line, "%s \"%s\" with %s", cmd, name, pstyle) ;
    } else {
        (void)sprintf(line, "%s \"%s\" title \"%s\" with %s", cmd, name,
                      title, pstyle) ;
    }
    
    /* send command to gnuplot  */
    gnuplot_cmd(line) ;
    nplots++ ;
    return ;
}



/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_plot1d_var2
 @memo		Plot a 2d graph from a list of dpoint.
 @param	d			Pointer to a list of doubles.
 @param	n_points	Number of doubles in the list.
 @param	title		Title of the plot.
 @return	void
 @doc
 
 Plots out a 2d graph from a list of dpoints. A dpoint is a struct
 containing two fields x and y (doubles) which are plotted as they
 are on the gnuplot session.
 
 \begin{verbatim}
 gnuplot_ctrl    *h ;
 dpoint          d[50] ;
 int             i ;
 
 h = gnuplot_init() ;
 for (i=0 ; i<50 ; i++) {
 d[i].x = (double)(i)/10.0 ;
 d[i].y = d[i].x * d[i].x ;
 }
 gnuplot_plot1d_var2(h, d, 50, "parabola") ;
 sleep(2) ;
 gnuplot_close(h) ;
 \end{verbatim}
 */
/*--------------------------------------------------------------------------*/

inline void Gnuplot_Control::gnuplot_plot1d_var2(
                                          dpoint          *   d,
                                          int                 n_points,
                                          char            *   title
)
{
    int         i ;
    FILE    *   tmp ;
    char    *   name ;
    char        cmd[GP_CMD_SIZE] ;
    char        line[GP_CMD_SIZE] ;
    
    /* can we open one more temporary file? */
    if (ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return ;
    }
    
    /* Open temporary file for output   */
    if ((name = tmpnam(NULL)) == (char*)NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }
    if ((tmp = fopen(name, "w")) == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }
    
    /* Store file name in array for future deletion */
    (void)strcpy(to_delete[ntmp], name) ;
    ntmp ++ ;
    
    /* Write data to this file  */
    for (i=0 ; i<n_points ; i++) {
        (void)fprintf(tmp, "%g %g\n", d[i].x, d[i].y) ;
    }
    (void)fflush(tmp) ;
    (void)fclose(tmp) ;

    
    /* Command to be sent to gnuplot    */
    if (nplots > 0) {
        (void)strcpy(cmd, "replot") ;
    } else {
        (void)strcpy(cmd, "plot") ;
    }
    
    if (title == NULL) {
        (void)sprintf(line, "%s \"%s\" with %s", cmd, name, pstyle) ;
    } else {
        (void)sprintf(line, "%s \"%s\" title \"%s\" with %s", cmd, name,
                      title, pstyle) ;
    }
    
    /* send command to gnuplot  */
    gnuplot_cmd(line) ;
    nplots++ ;

    return ;
}


inline void Gnuplot_Control::gnuplot_plot1d_var2(
                                                 double          *   x,
                                                 double          *   y,
                                                 int                 n_points,
                                                 char            *   title
)
{
    int         i ;
    FILE    *   tmp ;
    char    *   name ;
    char        cmd[GP_CMD_SIZE] ;
    char        line[GP_CMD_SIZE] ;
    
    /* can we open one more temporary file? */
    if (ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return ;
    }
    
    /* Open temporary file for output   */
    if ((name = tmpnam(NULL)) == (char*)NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }
    if ((tmp = fopen(name, "w")) == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }
    
    /* Store file name in array for future deletion */
    (void)strcpy(to_delete[ntmp], name) ;
    ntmp ++ ;
    
    /* Write data to this file  */
    for (i=0 ; i<n_points ; i++) {
        (void)fprintf(tmp, "%g %g\n", x[i], y[i]) ;
    }
    (void)fflush(tmp) ;
    (void)fclose(tmp) ;
    
    
    /* Command to be sent to gnuplot    */
    if (nplots > 0) {
        (void)strcpy(cmd, "replot") ;
    } else {
        (void)strcpy(cmd, "plot") ;
    }
    
    if (title == NULL) {
        (void)sprintf(line, "%s \"%s\" with %s", cmd, name, pstyle) ;
    } else {
        (void)sprintf(line, "%s \"%s\" title \"%s\" with %s", cmd, name,
                      title, pstyle) ;
    }
    
    /* send command to gnuplot  */
    gnuplot_cmd(line) ;
    nplots++ ;
    
    return ;
}





/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_plot_slope
 @memo		Plot a slope on a gnuplot session.
 @param	handle		Gnuplot session control handle.
 @param	a			Slope.
 @param	b			Intercept.
 @param	title		Title of the plot.
 @return	void
 @doc
 
 Plot a slope on a gnuplot session. The provided slope has an
 equation of the form:
 
 \begin{verbatim}
 y = ax+b
 \end{verbatim}
 
 Example:
 
 \begin{verbatim}
 gnuplot_ctrl        h ;
 double              a, b ;
 
 h.gnuplot_init() ;
 gnuplot_plot_slope(h, 1.0, 0.0, "unity slope") ;
 sleep(2) ;
 gnuplot_close(h) ;
 \end{verbatim}
 
 */
/*--------------------------------------------------------------------------*/


inline void Gnuplot_Control::gnuplot_plot_slope(
                                         double              a,
                                         double              b,
                                         char            *   title
)
{
    char    stitle[GP_TITLE_SIZE] ;
    char    cmd[GP_CMD_SIZE] ;
    
    if (title == NULL) {
        (void)strcpy(stitle, "no title") ;
    } else {
        (void)strcpy(stitle, title) ;
    }
    
    if (nplots > 0) {
        (void)sprintf(cmd, "replot %g * x + %g title \"%s\" with %s",
                      a, b, title, pstyle) ;
    } else {
        (void)sprintf(cmd, "plot %g * x + %g title \"%s\" with %s",
                      a, b, title, pstyle) ;
    }
    gnuplot_cmd(cmd) ;
    nplots++ ;
    return ;
}



inline void Gnuplot_Control::gnuplot_plot_equation (char *equation, char *title )
/*-------------------------------------------------------------------------*/
/**
 @name		gnuplot_plot_equation
 @memo		Plot a curve of given equation y=f(x).
 @param	equation	Equation to plot.
 @param	title		Title of the plot.
 @return	void
 @doc
 
 Plots out a curve of given equation. The general form of the
 equation is y=f(x), you only provide the f(x) side of the equation.
 
 Example:
 
 \begin{verbatim}
 gnuplot_ctrl    *h ;
 char            eq[80] ;
 
 h = gnuplot_init() ;
 strcpy(eq, "sin(x) * cos(2*x)") ;
 gnuplot_plot_equation(h, eq, "sine wave", normal) ;
 gnuplot_close(h) ;
 \end{verbatim}
 
 */
/*--------------------------------------------------------------------------*/

{
    char  cmd[GP_CMD_SIZE];
    char  plot_str[GP_EQ_SIZE];
    char  title_str[GP_TITLE_SIZE];
    
    if ( title == NULL )
    {
        ( void ) strcpy ( title_str, "no title" );
    } 
    else
    {
        ( void ) strcpy ( title_str, title);
    }
    
    if (nplots > 0 )
    {
        ( void ) strcpy ( plot_str, "replot" );
    } 
    else
    {
        ( void ) strcpy ( plot_str, "plot" );
    }
    
    ( void ) sprintf ( cmd, "%s %s title \"%s\" with %s", 
                      plot_str, equation, title_str, pstyle );
    
    gnuplot_cmd (cmd );
    
    nplots++;
    
    return;
}




#endif // _GNUPLOT_INCLUDED
