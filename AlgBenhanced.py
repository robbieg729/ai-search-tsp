############
############ ALTHOUGH I GIVE YOU THIS TEMPLATE PROGRAM WITH THE NAME 'skeleton.py', 
############ YOU CAN RENAME IT TO ANYTHING YOU LIKE. HOWEVER, FOR THE PURPOSES OF 
############ THE EXPLANATION IN THESE COMMENTS, I ASSUME THAT THIS PROGRAM IS STILL 
############ CALLED 'skeleton.py'.
############
############ IF YOU WISH TO IMPORT STANDARD MODULES, YOU CAN ADD THEM AFTER THOSE BELOW.
############ NOTE THAT YOU ARE NOT ALLOWED TO IMPORT ANY NON-STANDARD MODULES! TO SEE
############ THE STANDARD MODULES, TAKE A LOOK IN 'validate_before_handin.py'.
############

import os
import sys
import time
import random as rd
import copy

############
############ NOW PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS.
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############ BY 'DO NOT TOUCH' I REALLY MEAN THIS. EVEN CHANGING THE SYNTAX, BY
############ ADDING SPACES OR COMMENTS OR LINE RETURNS AND SO ON, COULD MEAN THAT
############ CODES WILL NOT RUN WHEN I RUN THEM!
############

def read_file_into_string(input_file, ord_range):
    the_file = open(input_file, 'r') 
    current_char = the_file.read(1) 
    file_string = ""
    length = len(ord_range)
    while current_char != "":
        i = 0
        while i < length:
            if ord(current_char) >= ord_range[i][0] and ord(current_char) <= ord_range[i][1]:
                file_string = file_string + current_char
                i = length
            else:
                i = i + 1
        current_char = the_file.read(1)
    the_file.close()
    return file_string

def remove_all_spaces(the_string):
    length = len(the_string)
    new_string = ""
    for i in range(length):
        if the_string[i] != " ":
            new_string = new_string + the_string[i]
    return new_string

def integerize(the_string):
    length = len(the_string)
    stripped_string = "0"
    for i in range(0, length):
        if ord(the_string[i]) >= 48 and ord(the_string[i]) <= 57:
            stripped_string = stripped_string + the_string[i]
    resulting_int = int(stripped_string)
    return resulting_int

def convert_to_list_of_int(the_string):
    list_of_integers = []
    location = 0
    finished = False
    while finished == False:
        found_comma = the_string.find(',', location)
        if found_comma == -1:
            finished = True
        else:
            list_of_integers.append(integerize(the_string[location:found_comma]))
            location = found_comma + 1
            if the_string[location:location + 5] == "NOTE=":
                finished = True
    return list_of_integers

def build_distance_matrix(num_cities, distances, city_format):
    dist_matrix = []
    i = 0
    if city_format == "full":
        for j in range(num_cities):
            row = []
            for k in range(0, num_cities):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    elif city_format == "upper_tri":
        for j in range(0, num_cities):
            row = []
            for k in range(j):
                row.append(0)
            for k in range(num_cities - j):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    else:
        for j in range(0, num_cities):
            row = []
            for k in range(j + 1):
                row.append(0)
            for k in range(0, num_cities - (j + 1)):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    if city_format == "upper_tri" or city_format == "strict_upper_tri":
        for i in range(0, num_cities):
            for j in range(0, num_cities):
                if i > j:
                    dist_matrix[i][j] = dist_matrix[j][i]
    return dist_matrix

def read_in_algorithm_codes_and_tariffs(alg_codes_file):
    flag = "good"
    code_dictionary = {}   
    tariff_dictionary = {}  
    if not os.path.exists(alg_codes_file):
        flag = "not_exist"  
        return code_dictionary, tariff_dictionary, flag
    ord_range = [[32, 126]]
    file_string = read_file_into_string(alg_codes_file, ord_range)  
    location = 0
    EOF = False
    list_of_items = []  
    while EOF == False: 
        found_comma = file_string.find(",", location)
        if found_comma == -1:
            EOF = True
            sandwich = file_string[location:]
        else:
            sandwich = file_string[location:found_comma]
            location = found_comma + 1
        list_of_items.append(sandwich)
    third_length = int(len(list_of_items)/3)
    for i in range(third_length):
        code_dictionary[list_of_items[3 * i]] = list_of_items[3 * i + 1]
        tariff_dictionary[list_of_items[3 * i]] = int(list_of_items[3 * i + 2])
    return code_dictionary, tariff_dictionary, flag

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY!
############
############ THE RESERVED VARIABLE 'input_file' IS THE CITY FILE UNDER CONSIDERATION.
############
############ IT CAN BE SUPPLIED BY SETTING THE VARIABLE BELOW OR VIA A COMMAND-LINE
############ EXECUTION OF THE FORM 'python skeleton.py city_file.txt'. WHEN SUPPLYING
############ THE CITY FILE VIA A COMMAND-LINE EXECUTION, ANY ASSIGNMENT OF THE VARIABLE
############ 'input_file' IN THE LINE BELOW iS SUPPRESSED.
############
############ IT IS ASSUMED THAT THIS PROGRAM 'skeleton.py' SITS IN A FOLDER THE NAME OF
############ WHICH IS YOUR USER-NAME, E.G., 'abcd12', WHICH IN TURN SITS IN ANOTHER
############ FOLDER. IN THIS OTHER FOLDER IS THE FOLDER 'city-files' AND NO MATTER HOW
############ THE NAME OF THE CITY FILE IS SUPPLIED TO THIS PROGRAM, IT IS ASSUMED THAT 
############ THE CITY FILE IS IN THE FOLDER 'city-files'.
############

input_file = "AISearchfile012.txt"

############
############ PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS.
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############

if len(sys.argv) > 1:
    input_file = sys.argv[1]

##### begin change 1 #####
the_particular_city_file_folder = "city-files"
path_for_city_files = "../" + the_particular_city_file_folder
##### end change 1   #####
    
if os.path.isfile(path_for_city_files + "/" + input_file):
    ord_range = [[32, 126]]
    file_string = read_file_into_string(path_for_city_files + "/" + input_file, ord_range)
    file_string = remove_all_spaces(file_string)
    print("I have found and read the input file " + input_file + ":")
else:
    print("*** error: The city file " + input_file + " does not exist in the folder '" + the_particular_city_file_folder + "'.")
    sys.exit()

location = file_string.find("SIZE=")
if location == -1:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()
    
comma = file_string.find(",", location)
if comma == -1:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()
    
num_cities_as_string = file_string[location + 5:comma]
num_cities = integerize(num_cities_as_string)
print("   the number of cities is stored in 'num_cities' and is " + str(num_cities))

comma = comma + 1
stripped_file_string = file_string[comma:]
distances = convert_to_list_of_int(stripped_file_string)

counted_distances = len(distances)
if counted_distances == num_cities * num_cities:
    city_format = "full"
elif counted_distances == (num_cities * (num_cities + 1))/2:
    city_format = "upper_tri"
elif counted_distances == (num_cities * (num_cities - 1))/2:
    city_format = "strict_upper_tri"
else:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()

dist_matrix = build_distance_matrix(num_cities, distances, city_format)
print("   the distance matrix 'dist_matrix' has been built.")

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY!
############
############ YOU NOW HAVE THE NUMBER OF CITIES STORED IN THE INTEGER VARIABLE 'num_cities'
############ AND THE TWO_DIMENSIONAL MATRIX 'dist_matrix' HOLDS THE INTEGER CITY-TO-CITY 
############ DISTANCES SO THAT 'dist_matrix[i][j]' IS THE DISTANCE FROM CITY 'i' TO CITY 'j'.
############ BOTH 'num_cities' AND 'dist_matrix' ARE RESERVED VARIABLES AND SHOULD FEED
############ INTO YOUR IMPLEMENTATIONS.
############

############
############ THERE NOW FOLLOWS CODE THAT READS THE ALGORITHM CODES AND TARIFFS FROM
############ THE TEXT-FILE 'alg_codes_and_tariffs.txt' INTO THE RESERVED DICTIONARIES
############ 'code_dictionary' AND 'tariff_dictionary'. DO NOT AMEND THIS CODE!
############ THE TEXT FILE 'alg_codes_and_tariffs.txt' SHOULD BE IN THE SAME FOLDER AS
############ THE FOLDER 'city-files' AND THE FOLDER WHOSE NAME IS YOUR USER-NAME, E.G., 'abcd12'.
############

##### begin change 2 #####
the_particular_alg_codes_and_tariffs = "alg_codes_and_tariffs.txt"
path_for_alg_codes_and_tariffs = "../" + the_particular_alg_codes_and_tariffs
##### end change 2   #####

code_dictionary, tariff_dictionary, flag = read_in_algorithm_codes_and_tariffs(path_for_alg_codes_and_tariffs)

if flag != "good":
    print("*** error: The text file 'alg_codes_and_tariffs.txt' does not exist.")
    sys.exit()

print("The codes and tariffs have been read from 'alg_codes_and_tariffs.txt':")

############
############ HAVE YOU TOUCHED ANYTHING ABOVE? BECAUSE EVEN CHANGING ONE CHARACTER OR
############ ADDING ONE SPACE OR LINE RETURN WILL MEAN THAT THE PROGRAM YOU HAND IN
############ MIGHT NOT RUN PROPERLY! SORRY TO GO ON ABOUT THIS BUT YOU NEED TO BE 
############ AWARE OF THIS FACT!
############
############ YOU NOW NEED TO SUPPLY SOME PARAMETERS.
############
############ THE RESERVED STRING VARIABLE 'my_user_name' SHOULD BE SET AT YOUR
############ USER-NAME, E.G., "abcd12"
############

my_user_name = "cnkz75"

############
############ YOU CAN SUPPLY, IF YOU WANT, YOUR FULL NAME. THIS IS NOT USED AT ALL BUT SERVES AS
############ AN EXTRA CHECK THAT THIS FILE BELONGS TO YOU. IF YOU DO NOT WANT TO SUPPLY YOUR
############ NAME THEN EITHER SET THE STRING VARIABLES 'my_first_name' AND 'my_last_name' AT 
############ SOMETHING LIKE "Mickey" AND "Mouse" OR AS THE EMPTY STRING (AS THEY ARE NOW;
############ BUT PLEASE ENSURE THAT THE RESERVED VARIABLES 'my_first_name' AND 'my_last_name'
############ ARE SET AT SOMETHING).
############

my_first_name = "Robbie"
my_last_name = "Goodall"

############
############ YOU NEED TO SUPPLY THE ALGORITHM CODE IN THE RESERVED STRING VARIABLE 'algorithm_code'
############ FOR THE ALGORITHM YOU ARE IMPLEMENTING. IT NEEDS TO BE A LEGAL CODE FROM THE TEXT-FILE
############ 'alg_codes_and_tariffs.txt' (READ THIS FILE TO SEE THE CODES).
############

algorithm_code = "AC"

############
############ DO NOT TOUCH OR ALTER THE CODE BELOW! YOU HAVE BEEN WARNED!
############

if not algorithm_code in code_dictionary:
    print("*** error: the algorithm code " + algorithm_code + " is illegal")
    sys.exit()
print("   your algorithm code is legal and is " + algorithm_code + " -" + code_dictionary[algorithm_code] + ".")

############
############ YOU CAN ADD A NOTE THAT WILL BE ADDED AT THE END OF THE RESULTING TOUR FILE IF YOU LIKE,
############ E.G., "in my basic greedy search, I broke ties by always visiting the first 
############ city found" BY USING THE RESERVED STRING VARIABLE 'added_note' OR LEAVE IT EMPTY
############ IF YOU WISH. THIS HAS NO EFFECT ON MARKS BUT HELPS YOU TO REMEMBER THINGS ABOUT
############ YOUR TOUR THAT YOU MIGHT BE INTERESTED IN LATER.
############

added_note = ""

############
############ NOW YOUR CODE SHOULD BEGIN.
############

start = time.time()

def compute_nearest_neighbor_tour(dist_matrix, start_city):

    # initialize the tour, tour length, and list of unvisited cities
    tour = [start_city]
    tour_length = 0
    unvisited_cities = list(range(num_cities))
    unvisited_cities.remove(start_city)
    
    # create a loop which will execute num_cities - 1 times
    for i in range(0, num_cities - 1):

        # get current city
        current_city = tour[-1]

        # get the row of the distance matrix corresponding to the current city
        row = dist_matrix[current_city]

        # get the distances from the current city to all unvisited cities
        legal_distances = [row[city] for city in unvisited_cities]

        # get next city by getting the index of the minimum value along the row of the current city, only checking the columns 
        # for the unvisited cities. Then use that index to find the corresponding city
        next_city = unvisited_cities[legal_distances.index(min(legal_distances))]

        # add the next city to the tour, and remove it from the unvisited cities list
        tour.append(next_city)
        unvisited_cities.remove(next_city)

        # add the distance from the current city to the next city, to the tour length
        tour_length += dist_matrix[current_city][next_city]

        # if this is the final iteration of the loop, also add the distance from the end city to the start city
        if i == num_cities - 2:
            tour_length += dist_matrix[next_city][start_city]

    # return the tour and the tour length
    return [tour, tour_length]

# get the tour and tour length of the nearest neighbor tour, starting from 0, and set them to be the current best
tour, tour_length = compute_nearest_neighbor_tour(dist_matrix, 0)

# set parameter values alpha, beta, rho, N, max_iter, p_best, p_dec, tao_max, and tao_min. Note p_best and p_dec are only used
# to calculate the value of tao_min
alpha = 1
beta = 2
rho = 0.3
N = num_cities
max_iter = 10000
p_best = 5e-3
p_dec = (p_best) ** (1 / num_cities)
tao_max = (1 / rho) * (1 / tour_length)
tao_min = (tao_max * (1 - p_dec)) / ((num_cities / 2 - 1) * p_dec)

# class to store information about a single ant
class Ant:

    def __init__(self, start_city):

        # set starting city
        self.start_city = start_city

        # set all variables to their base values (what they should be at the start of an iteration)
        self.reset_variables()

    def select_next_city(self):

        # if there is only one unvisited city left, return that city
        if len(self.unvisited_cities) == 1:
            return self.unvisited_cities[0]

        # else, select the next city using probabilistic methods
        else:

            # most recently visited city in the ant's tour
            current_city = self.tour[-1]

            # get the row of the distance matrix corresponding to the current city and change all values of 0 to 0.00001
            # to avoid division by zero errors
            row = [1e-5 if d == 0 else d for d in dist_matrix[current_city]]

            # get 'probabilities' of all unvisited cities being the next selected city. Note they are not actual probabilities
            # they are just weights which can be used as probabilities
            probabilities = [(pheromone_levels[current_city][city] ** alpha) * ((1 / row[city]) ** beta) for city in self.unvisited_cities]

            # select and return the next city based on a random sample from the unvisited cities, and their probabilities
            next_city = rd.choices(self.unvisited_cities, weights=probabilities, k=1)[0]
            return next_city

    def reset_variables(self):
        # reset necessary variables after an iteration or during initialization

        # set list of unvisited cities to be all cities, with the starting city removed
        self.unvisited_cities = list(range(num_cities))
        self.unvisited_cities.remove(self.start_city)

        # set tour length to be 0 and current tour to contain only the starting city
        self.tour_len = 0
        self.tour = [self.start_city]

    def copy_info_from_ant(self, ant):
        # copy another ant's info to this ant
        self.tour = ant.tour
        self.tour_len = ant.tour_len
        self.start_city = ant.start_city

# initialize pheromone levels at each edge to be tao_zero
pheromone_levels = [[tao_max for i in range(num_cities)] for j in range(num_cities)]

# get initial list of ants
ants = [Ant(city) for city in range(N)]

# initialize the global best ant to store the information from the tour produced by the nearest neighbors algorithm
global_best_ant = Ant(0)
global_best_ant.tour = tour
global_best_ant.tour_len = tour_length
global_best_ant.unvisited_cities = []

# variable which when true will end the computations early for time constraints
end_iterations = False 

for i in range(0, max_iter):

    # initialize the ant with the best tour from this iteration to be an Ant with an incomplete tour
    iteration_best_ant = Ant(0)

    # for each ant
    for ant in ants:

        # if there is not enough time to keep simulating, break the loop and set end_iterations to True
        if time.time() - start > 59:
            end_iterations = True
            break

        # while the ant's tour is not complete
        while len(ant.unvisited_cities) > 0:

            # select the next city based on probabilistic measures
            next_city = ant.select_next_city()

            # add the length of the chosen edge to the ant's tour length
            ant.tour_len += int(dist_matrix[ant.tour[-1]][next_city])

            # add the chosen city to the ant's tour
            ant.tour.append(next_city)

            # remove the chosen city from the ant's unvisited cities list
            ant.unvisited_cities.remove(next_city)
            
        # once the tour is complete, add the distance from the final city to the start city, to the ant's tour length
        ant.tour_len += int(dist_matrix[ant.tour[-1]][ant.start_city])
        
        # if the tour length of the ant is shorter than the current best, update the tour and tour length for the current best.
        # Also update the information for the global best ant, and update the value for tao_max since a new best tour has been
        # found
        if ant.tour_len < tour_length:
            tour_length = ant.tour_len
            tour = ant.tour
            global_best_ant.copy_info_from_ant(ant)
            tao_max = (1 / rho) * (1 / tour_length)
            
        # if the tour length of the ant is shorter than the current best from the current iteration, or this is the first ant
        # in the iteration, update the information for the iteration best ant
        if len(iteration_best_ant.tour) != num_cities or ant.tour_len < iteration_best_ant.tour_len:
            iteration_best_ant.copy_info_from_ant(ant)

        # reset the ant's variables, such as tour, tour_len, and unvisited_cities, in preparation for the next iteration
        ant.reset_variables()

    # break the loop if we need to terminate the program early due to time constraints
    if end_iterations == True:
        break
    
    # evaporate pheromone from all edges, capping at tao_min if necessary
    new_pheromone_levels = [[(1 - rho) * p if (1 - rho) * p > tao_min else tao_min for p in pheromone_levels[j]] for j in range(num_cities)]

    # set the ant which will add pheromone to the iteration best ant. For the first 25 iterations we use this ant. After that,
    # we use the global best ant according to the rules below
    pheromone_updating_ant = iteration_best_ant
    
    # After 25 iterations, we use the global best ant at different frequencies. For iterations 25-75, we use it every 5
    # iterations. For iterations 76-125, we use it every 3. For iterations 126-250, we use it every 2, and after that point,
    # we use the global best ant in every iteration
    if i >= 25:
        global_best_freq = 5
        if i >= 75 and i < 125:
            global_best_freq = 3
        elif i >= 125 and i < 250:
            global_best_freq = 2
        else:
            global_best_freq = 1
        if i % global_best_freq == 0:
            pheromone_updating_ant = global_best_ant

    # let the ant deposit pheromone onto its trail, making sure to cap all trails at tao_max
    for j in range(0, num_cities):
        current_city = pheromone_updating_ant.tour[j]
        next_city = pheromone_updating_ant.tour[0] if j == num_cities - 1 else pheromone_updating_ant.tour[j + 1]
        new_pheromone_level = tao_max if new_pheromone_levels[current_city][next_city] + 1 / pheromone_updating_ant.tour_len > tao_max else new_pheromone_levels[current_city][next_city] + 1 / pheromone_updating_ant.tour_len
        new_pheromone_levels[current_city][next_city] = new_pheromone_level

    # copy the new pheromone levels to the current pheromone levels
    pheromone_levels = copy.deepcopy(new_pheromone_levels)

# ensure the tour and tour lengths are integers. This shouldn't change the actual values of these variables, only their 
# types, as the code above may store them as floats (or a list of floats). We need them to be integers because of the check
# in the code later on 
tour = [int(city) for city in tour]
tour_length = int(tour_length)

















############
############ YOUR CODE SHOULD NOW BE COMPLETE AND WHEN EXECUTION OF THIS PROGRAM 'skeleton.py'
############ REACHES THIS POINT, YOU SHOULD HAVE COMPUTED A TOUR IN THE RESERVED LIST VARIABLE 'tour', 
############ WHICH HOLDS A LIST OF THE INTEGERS FROM {0, 1, ..., 'num_cities' - 1} SO THAT EVERY INTEGER
############ APPEARS EXACTLY ONCE, AND YOU SHOULD ALSO HOLD THE LENGTH OF THIS TOUR IN THE RESERVED
############ INTEGER VARIABLE 'tour_length'.
############

############
############ YOUR TOUR WILL BE PACKAGED IN A TOUR FILE OF THE APPROPRIATE FORMAT AND THIS TOUR FILE'S,
############ NAME WILL BE A MIX OF THE NAME OF THE CITY FILE, THE NAME OF THIS PROGRAM AND THE
############ CURRENT DATA AND TIME. SO, EVERY SUCCESSFUL EXECUTION GIVES A TOUR FILE WITH A UNIQUE
############ NAME AND YOU CAN RENAME THE ONES YOU WANT TO KEEP LATER.
############

############
############ DO NOT TOUCH OR ALTER THE CODE BELOW THIS POINT! YOU HAVE BEEN WARNED!
############

flag = "good"
length = len(tour)
for i in range(0, length):
    if isinstance(tour[i], int) == False:
        flag = "bad"
    else:
        tour[i] = int(tour[i])
if flag == "bad":
    print("*** error: Your tour contains non-integer values.")
    sys.exit()
if isinstance(tour_length, int) == False:
    print("*** error: The tour-length is a non-integer value.")
    sys.exit()
tour_length = int(tour_length)
if len(tour) != num_cities:
    print("*** error: The tour does not consist of " + str(num_cities) + " cities as there are, in fact, " + str(len(tour)) + ".")
    sys.exit()
flag = "good"
for i in range(0, num_cities):
    if not i in tour:
        flag = "bad"
if flag == "bad":
    print("*** error: Your tour has illegal or repeated city names.")
    sys.exit()
check_tour_length = 0
for i in range(0, num_cities - 1):
    check_tour_length = check_tour_length + dist_matrix[tour[i]][tour[i + 1]]
check_tour_length = check_tour_length + dist_matrix[tour[num_cities - 1]][tour[0]]
if tour_length != check_tour_length:
    flag = print("*** error: The length of your tour is not " + str(tour_length) + "; it is actually " + str(check_tour_length) + ".")
    sys.exit()
print("You, user " + my_user_name + ", have successfully built a tour of length " + str(tour_length) + "!")

local_time = time.asctime(time.localtime(time.time()))
output_file_time = local_time[4:7] + local_time[8:10] + local_time[11:13] + local_time[14:16] + local_time[17:19]
output_file_time = output_file_time.replace(" ", "0")
script_name = os.path.basename(sys.argv[0])
if len(sys.argv) > 2:
    output_file_time = sys.argv[2]
output_file_name = script_name[0:len(script_name) - 3] + "_" + input_file[0:len(input_file) - 4] + "_" + output_file_time + ".txt"

f = open(output_file_name,'w')
f.write("USER = " + my_user_name + " (" + my_first_name + " " + my_last_name + "),\n")
f.write("ALGORITHM CODE = " + algorithm_code + ", NAME OF CITY-FILE = " + input_file + ",\n")
f.write("SIZE = " + str(num_cities) + ", TOUR LENGTH = " + str(tour_length) + ",\n")
f.write(str(tour[0]))
for i in range(1,num_cities):
    f.write("," + str(tour[i]))
f.write(",\nNOTE = " + added_note)
f.close()
print("I have successfully written your tour to the tour file:\n   " + output_file_name + ".")
    
    











    


