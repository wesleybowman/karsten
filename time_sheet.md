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

```
