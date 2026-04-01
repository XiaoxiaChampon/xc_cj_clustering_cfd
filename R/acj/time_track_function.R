############################################################
# Copyright 2023 Xiaoxia Champon

# Permission is hereby granted, free of charge, to any person 
# obtaining a copy of this software and associated documentation 
# files (the “Software”), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, 
# publish, distribute, sublicense, and/or sell copies of the Software, 
# and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
######################################################################
#
# Purpose: Functions to perform Categorical Functional Data Hypothesis Testing
#           File that contains all the functions necessary to generate data 
# Author:  Xiaoxia Champon
# Date: 10/26/2023
#
##############################################################


time_elapsed <<- list()
last_time <<- 0
row_name <<- NULL


#' Function to start recording the time for one task
#' @param task_name : task name that needs to track the time
timeKeeperStart <- function(task_name)
{
    row_name <<- task_name
    if(FALSE == row_name %in% names(time_elapsed))
    {
        time_elapsed[[row_name]] <<- NULL
    }
    last_time <<- Sys.time()
}


#' Function to print the time taken by the task 
timeKeeperNext <- function()
{
    this_time <- Sys.time()
    this_section_time <- this_time - last_time
    cat("\n--------------------\n",
        row_name, "\n\ttook:", capture.output(this_section_time), 
        "\n====================\n")
    time_elapsed[[row_name]] <<- append(time_elapsed[[row_name]], this_section_time)
    last_time <<- this_time
}
