```diff
- 2014/04/25
+ 09:00 - 12:00 -- 
    Got my computer to boot
+ 12:00 - 1:30 -- 
    Debugged the year long automatic script (while it was running), 
    and fixed the code on Github.
+ 2-3 --
    Went to lunch
+ 3-4 --
    Filled out payform and finished the booting stuff from the morning.
+ 4-5 --
    Looking at Plausible_Turbine_Locations.m, also asked Duane about
    software on the Acadia cluster, hoping to get git there.

    Began breakdown analysis of Plausible_turbine_locations, seeing what
    function calls are made, and going into those and seeing what is called
    until we get down to matlab function calls (except uTides).

    Most seem reasonable, need to figure out about uTides.
+5-5:30 --
    Worked a little on oct2py, and looked into IPython notebooks.

- 2014/04/27
+ 9:42-13:09 --
    Started conversion process on turbine array code.

    From doing the loading, seems like matlab takes 4.851s while python takes
    31.7ms. More testing may need to be done.

    Had some errors, but that was because MatLab starts at index 1 while python
    starts at index zero. The error came in because trinodes is an array of
    indices, so I had to subtract one from all of them for it to be correct.

- 2014/04/28
+ 8:30-9 --
    Looked at the year long automation, there is a bug I was trying to find,
    need Aidan to point it out. (Turned out there wasn't a bug, just the code
    Aidan was using was the updated version.)

+ 9-15:47 --
    More turbine array coding. Converted most of the code. Jon and Andy came
    by, so some discussion with them about linux and python.


    Still trying to get cf_u_rated_turbs to work. Having some issues with
    structures from Matlab. Got them to work, but now it has stopped working
    for some reason.

    calculate_power is the issue. Having difficulties seeing the problem since
    its all functional based. Wrote  test script to make it easier. Seems to be
    a broadcasting error with the structured numpy arrays. Will look more into
    it tomorrow.

- 2014/04/29
+ 9-10:31 --
    Ironed out the bugs in yesterdays code. Used recarrays instead of
    structured ones. Running both matlab and the code now to see if the returns
    are the same. No impressive speedups yet. If anything Matlab may be faster,
    but I see areas for big improvement.

+ 10:31-11:40 --
    Saved matlab turbines file so that I can review it without having to run
    it. Checked to see if code was doing correct calculations. Got rid of a for
    loop. Waiting to see if code is still correct.

+ 11:40-12:15 --
    Lunch

+ 12:30 - 13:30 --
    Still debugging. Let run and went to see Dr. K.

+ 13:30 -15:00 --
    Working on turning recarrays into pandas DataFrames (this is easy, one line
    command). Also trying to save as h5 (also easy) and then load those h5
    files into matlab (proving difficult).

    Converted turbines into DataFrame, then saved the dataframe as a csv, which
    I could then use to load the variables into matlab.

    Tidied up code.

+ 15-16:24 --
    Started generalRunFiles. Talked with Aidan about the things that need to be
    able to be changed.

    Got it up and working. Asks for 17 different variables, then creates the
    .nml and the .sh files needed for qsubbing.


```
