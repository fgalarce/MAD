#ifndef SUBWINDOW
#define SUBWINDOW


//#include <string>

/* <<show_subwindow()>> generates a subwindow within the viewport
 * By default shows the par editor. In the future (probably), will support other stuff*/
void show_subwindow(int win_type=0);

void display_ParEditor();

void display_GraphTest();

void display_DemoPlot();

#endif