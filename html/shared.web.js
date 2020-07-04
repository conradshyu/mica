var progressEnd         = 10;       // set to number of progress <span>'s.
var progressColor       = 'blue';   // set to progress bar color
var progressInterval    = 1000;     // set to time between updates (milli-seconds)
var progressAt          = 0;
var progressTimer;

function setPeriod( query_period )
{
    progressInterval = query_period;
}

function progressClear()
{
    for ( var i = 0; i < progressEnd; i++ )
    {
        document.getElementById( 'progress' + i ).style.backgroundColor =
            'transparent';
    }

    progressAt = 0;
}

function progressStop()
{
    clearTimeout( progressTimer );
    progressClear();
}

function progressUpdate()
{
    if ( progressAt < progressEnd )
    {
        document.getElementById( 'progress' + progressAt ).style.backgroundColor =
            progressColor;
        progressTimer = setTimeout( 'progressUpdate()', progressInterval );
        progressAt++;
    }
    else
    {
        progressStop();
    }
}

function popUp( file_name, width, height )
{
    day = new Date();
    id = day.getTime();

    win_parameter = 'toolbar=0, scrollbars=1, location=0, statusbar=1, menubar=0, resizable=1, width=' + width + ', height=' + height + ', left=200, top=200';
    eval( "page" + id + " = window.open( file_name, '" + id + "', '" + win_parameter + "');" );
}
