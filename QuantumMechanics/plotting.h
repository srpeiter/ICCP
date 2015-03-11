#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<random>

#ifndef PLOTTING_H
#define PLOTTING_H


class mydata
{
int particles;
double *xdata=NULL;
double *ydata=NULL;
double *zdata=NULL;
public:
mydata (double *xdata, int particles) : xdata(xdata), particles(particles) {}

mydata (double *xdata, double *ydata, int particles): xdata(xdata), ydata(ydata), particles(particles) {}

mydata (double *xdata, double *ydata, double *zdata, int particles): xdata(xdata), ydata(ydata), zdata(zdata), particles(particles) {}


void printtofile(char filename[]);
void plot2d();
void plot3d();
void histplot();
};

#endif

