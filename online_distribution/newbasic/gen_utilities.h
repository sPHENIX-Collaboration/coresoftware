/*
 *
 *-------------------------------------------------------------
 * general purpose utilities
 *-------------------------------------------------------------
 *
 *
 * This file contains some useful general purpose routines
 * 
 * haveCommand: reads a command from a stream 
 *              and returns an int indicating 
 *              whether there is a command waiting
 *
 */


extern int haveCommand( FILE *stream = stdin );
