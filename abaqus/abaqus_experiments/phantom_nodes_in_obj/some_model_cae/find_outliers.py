


with open("/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/phantom_nodes_in_obj_cae/some_model_cae/test.obj") as f:
    cur_line = f.readline() 
    line_num = 1
    while (cur_line != ""):
        line_components = cur_line.split(" ")  
        if line_components[0] == "v":
            # Check for outlier vertices.
            printed_line_num = False
            if (float(line_components[1]) < -1 or float(line_components[1]) > 16):
                print("Found on line number " + str(line_num) + ":")
                print("Unusual x value of " + str(line_components[1]) + ".")
                printed_line_num = True

            if (float(line_components[2]) < -1 or float(line_components[2]) > 16):
                if not printed_line_num:
                    print("Found on line number " + str(line_num) + ":")
                    printed_line_num = True
                print("Unusual y value of " + str(line_components[2]) + ".")

            if (float(line_components[3]) < -1 or float(line_components[3]) > 21):
                if not printed_line_num:
                    print("Found on line number " + str(line_num) + ":")
                printed_line_num = True
                print("Unusual z value of " + str(line_components[3]) + ".")
        cur_line = f.readline()
        line_num += 1
