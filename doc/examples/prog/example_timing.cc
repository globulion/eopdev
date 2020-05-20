clock_t t_time = -clock(); // initialize the clock
t_time += clock() // temporarily stop the clock in this place
t_time -= clock() // start again the clock here
cout << " o TIME LABEL: " << ((double)t_time*1000/CLOCKS_PER_SEC) << endl;
