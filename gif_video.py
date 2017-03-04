import os

def make_gif():
    directory = raw_input("Directory: ")
    dir_path = os.getcwd() + "/" + directory
    if os.path.exists(dir_path):
        output = dir_path.split("/")[-1]
        command = "convert -loop 0 {0}/*.png {0}/{1}.gif".format(dir_path,output)
        print command
        os.system(command)
    else:
        print "%s does not exist" % dir_path
    return

def make_video():
    directory = raw_input("Directory: ")
    dir_path = os.getcwd() + "/" + directory
    if os.path.exists(dir_path):
        output = dir_path.split("/")[-1]
        command = "ffmpeg -framerate 10 -pattern_type glob -i '{0}/*.png' {0}/{1}.mp4".format(dir_path,output)
        print command
        os.system(command)
    else:
        print "%s does not exist" % dir_path
    return

def main():
    n = 1
    while n == 1:
        prompt = raw_input("Gif/Video? (g/v): ")
        if prompt == "g":
            make_gif()
            n = 0
        if prompt == "v":
            make_video()
            n = 0
        elif prompt != "g" and prompt != "v":
            print "try again"

main()
