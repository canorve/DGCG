import subprocess

try:
    # pipe output to /dev/null for silence
    null = open("/dev/null", "w")
    subprocess.Popen("sextractor", stdout=null, stderr=null)
    null.close()

except OSError:
    print("git not found")
	sys.exit()

