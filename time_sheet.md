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

+ 16:24- 17:16 --
    Tested the general scripts. Had to wait on Aidan to check them. Still need
    to change bottom_roughness_minimum. Should also add making the RUN_PAR2 to
    the yearLongRuns as well.

- 2014/04/30
+ 8:56-9:52 --
    Working on bash line history in scripts. 

    Stopped since I was getting no where, and asked the internet for help.

+ 9:52-10:30 --
    Looking into timing python and matlab scripts that were written yesterday.
    Also, while those are running, testing more ways to handle structures from
    matlab. Specifically:
    - [ ] structured arrays
    - [ ] recarrays
    - [ ] pandas panels/dataframe
    - [ ] dictionaries

    Time for matlab is 1052s for cf_u_rated_turbs.
    Time for python is 2m52s.

+ 11:00-15:56 --
    Meeting with Dr. K, then lunch, then a second meeting. Talked about turbine
    analysis, and then validation of runs, and then UTides.

+ 15:56 - 16:30ish --
    Looking into U Tides.
    Running it on matlab to get coef.

    Computer crashed since I used all ram up.

- 2014/05/01
+ 7:45-9 --
    Run matlba code to get coef.

+ 9:10:46 --
    Trying to save pythonTurbines as an nc file and test to see that it loads
    into MatLab. Have Matlab coef and turbines saved to.mat files, and have
    python turbines saved so no code has to be run.

    Saving to nc requires all of the ram. I had 200 Mb when I started to save,
    and it jumped to 7.4 Gb to save. Crashed computer again, then redid it and
    got the file saved.
    When trying to load in the nc file, Matlab throws OutofMemory Error.

    Started workig on converting ut_solv.
    Got all converted up to ut_solv, and have done testing on some variables in
    ut_solv for comparison.

+ 10:46-13:00 --
    Starting UTide conversion (brace yourself).
    Got:
    - [x] ut_solv
    - [x] ut_solv1
    - [x] ut_slvinit

+ 13:00 - 13:27 --
    Lunch

+ 13:28-16:05 --
    Back to working on code, specifically ut_solv1.

    got ut_astron and am working on ut_constitsel. Problem right now is loading
    in the mat file. Wasn't the problem, multiplication was the problem. Used
    np.dot and it fixed the issue. Done for the day. Been looking at this for
    way too long today. Should start to go faster, since I am learning where
    the problems are.

- 2014/05/02
+ 14:00-16:00 --
    Conference call with Joel and the crew talking about GIS

- 2014/05/03
+ 11:00-11:23 --
    Emailed Codiga about UTides, and emailed Travis Oliphant about robustfit.
    I am pretty sure we don't need robustfit immediately, but for eventual
    completeness it will be nice.

+ 11:24-12:35 --
    Started converting code.
    Got to ut_FUV, and its a little longer, so lunch time.

+ 12:35 - 13:40 -- 
    Lunch

+ 13:40 -16:10 --
    Back to code.
    Got to the bottom on ut_FUV. Had problems with getting numbers to look the
    same. Usually small mistakes with the variables.

+ 16:10 - 16:19 --
    Got the issue. Indices can be a pain in the ass. Finished ut_FUV... not
    finished, but have F and U calculated. Still need to do V.

+ 16:20 -17:41 --
    Finishing up FUV. Finsihed. Worked when I did everything by hand, so once I
    run it all it should work.
    Got B finished, which is where we need to be. Will continue tomorrow.

- 2014/05/04
+ 9:36 -10:45 --
    Started working on code.

+ 14:00-14:30 --
    Code now produces correct coef.Lsmaj. Right now its not presented in the
    way MatLab does, but it does work.

+ 14:30 - 15:17 --
    Finished all Plausible_turbine_locations code. Now can be run in python.
    Still very rough, but is the backbone for more features to be added.
    Still need to get plots in.

- 2014/05/05
+ 08:15 - 15:48 --
    Started running the matlab plausible_turbine_locations to get time. Then am
    going to run the python version and compare. Also working on getting plots
    up.

    Did some plotting work, and got all things timed.
    Long meeting with Karston.
    Also worked on ut_reconstr and got a barebones version in.

- 2014/05/06
+ 9:18 -10:21 --
    Looked into some of Brian Polagye's code. Started writing the matlab/python
    comparison script.

+ 10:21 - 11:13 --
    Running matlab to get a coef to compare to.
    Got comparison to work.

    Editted the code to only show things that are not equal, that way its
    easier to see what is wrong.

+ Rest of day:
    Worked a bit on plotting the turbines, and trying to get them to work so we
    can make a movie.

- 2014/05/07
+ 08:30-09:37 --
    Starting making UTide github page. Made file hierarchy of how I want the
    package structured, this is tenative. Started __init__.py

+ 9:37-10:14 -- 
    Got the package structure going. Put the functions in their respective
    directories as specified by the hierarchy.

+ 10:14 - 12:24 --
    Packaging didn't work. Resorted back to old. made sure 1D data worked.
    Lunch.

+ 12:24 -15:09 --
    Confidence interval stuff.
    Got everything coded, just need to test all of it and get it working.

    Right now, pwelch is the problem.

+ 15:09-15:28 --
    Got welch to maybe work, havent compared. Waiting for code to run before
    doing comparison.

    Figured out a good way to debug. Much better than before. Much faster.

- 2014/05/08
+ 8:50 -9:13 --
    Making package work. Fast and easy solution, nothing complicated. K.I.S.S.

    from UTide import * now works with both ut_solv and ut_reconstr. There are
    some things to look into for the future, but for now it works as a package.
    It is not the most up to date version, but it is a working version for
    plausible_turbine_locations. May just start an unstable branch in that
    repository.

+ 9:13 -17:12  --
    Working on debugging the CI stuff.
    Got some of it done. Pwelch causing an issue. Worked for one of them, but
    for the second its off, and because its off, that error propagates, since
    we use that value in means and whatnot.

    Had a 2h K meeting with Joel and the crew. 

    Looked some more at CI. Can't see what is wrong.

+ 17:13 - 17:58 --
    Looked some more into CI, and then figured out what was causing differences
    between frq, turns out matlab reads in 4 decimals and python reads in more
    than that, causing the small discrepancy.

+ GRADUATION, missed a few days

+ 2014/05/15
- 9:00 - 9:23 --
    Got ADCP data and signed the doodle for the volunteering at the conference.

- 9:23 - --


```
