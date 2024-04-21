def parse_file(fname):
    fileraw = open(fname, 'r').read()

    filestr = uncomment_file(fileraw)

    data = {}

    i = 0
    sp = filestr.split(';')

    while i < len(sp):
        line = sp[i].replace(' ', '')
        if '{' in line:
            key = line.split('={')[0]
            newline = '={'.join(line.split('={')[1:])
            data[key], i, sp = read_level(newline, i, sp)
        elif line != '':
            raise RuntimeError('Unable to parse: %s' % line)

        i = i + 1

    return data


def read_level(line, i, sp):
    subdata = {}

    while True:
        if '{' in line:
            key = line.split('={')[0]
            newline = '={'.join(line.split('={')[1:])
            subdata[key], i, sp = read_level(newline, i, sp)
        elif '}' in line:
            return subdata, i, sp
        elif '=' in line:
            [key, value] = line.split('=')
            subdata[key] = value
        elif line != '':
            raise RuntimeError('Unable to parse: %s' % line)

        i = i + 1

        if i == len(sp):
            raise RuntimeError('EOF reached while parsing an input block. Did you forget to close a bracket?')

        line = sp[i].replace(' ', '')


# parse list, if length is provided the behavior is as follows:
# if parsed_length == 1, a list with repeating value is created
# if parsed_length > 1, the length of the list is checked

def parse_list(lst, typ=str, length=None):
    if type(lst) is str:
        lst = lst.strip('[').strip(']').replace(' ', '').split(',')

    if length is not None:
        if len(lst) == 1:
            lst = [lst[0] for i in range(length)]
        assert len(lst)==length, f"expected list of length 1 or {length}, got {len(lst)}" 

    return list(map(typ, lst))


def uncomment_file(fileraw):
    filestr = ''
    comment_mode = False

    # Go through all lines in the raw file
    for line in fileraw.split('\n'):

        # Check if we are in comment mode
        if comment_mode:

            # If so, try to find '*/' and remove only the part before it
            end = line.find('*/')
            if end >= 0:
                comment_mode = False
                line = line[end + len('*/'):]

            # If comment_mode is still enabled, don't include anything
            else:
                line = ''

        if not comment_mode:

            # If we are not in comment mode, remove all full comments from the line
            line = uncomment_line(line)

            # If there is a '/*' in the line, enable comment mode, and remove the part after '/*' from the line
            start = line.find('/*')
            if start >= 0:
                comment_mode = True
                line = line[:start]

        # Add the line to the file string
        filestr += line.replace('\t', '')

    return filestr


def uncomment_line(line):
    clean_line = line

    # Remove all comments from the line (assuming that comment_mode is False)
    while True:
        # If the first identifier is a '//', remove everything after it
        start_oneline = clean_line.find('//')
        start_block = clean_line.find('/*')
        if start_oneline >= 0:
            if start_oneline < start_block or start_block < 0:
                clean_line = clean_line[:start_oneline]

        # Remove everything in the first one-line block comment that is found
        start = clean_line.find('/*')
        end = clean_line.find('*/')
        if 0 <= start < end:
            clean_line = clean_line[:start] + clean_line[end + len('*/'):]
        else:
            # Exit if no comments are left
            break

    return clean_line


def soft_cast(value, typ):
    # This function attempts to convert value to typ
    # If this conversion fails, it returns the original value
    try:
        return typ(value)
    except:
        return value


def evaluate(value, coords, rank, extra_dict=None):
    # This function does a string evaluation of value, if possible
    if type(value) is str:
        eval_dict = get_eval_dict(coords, rank, extra_dict)
        return eval(value, {}, eval_dict)
    else:
        return value


def get_eval_dict(coords, rank, extra_dict=None):
    # This function builds a dictionary with the x,y,z coordinates of coords
    eval_dict = {'x': coords[0]}

    if rank >= 2:
        eval_dict.update({'y': coords[1]})
    if rank == 3:
        eval_dict.update({'z': coords[2]})

    # Add the extra dictionary if applicable
    if extra_dict is not None:
        eval_dict.update(extra_dict)

    return eval_dict
