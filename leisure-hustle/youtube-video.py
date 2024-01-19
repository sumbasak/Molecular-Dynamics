from yt_dlp import YoutubeDL

playlist_url = "https://www.youtube.com/playlist?list=PLqj1RDXdueExLtODAEs-3OWDql84EZjWi"

ydl_opts = {
    'format': 'bestaudio/best',
    'postprocessors': [{
        'key': 'FFmpegExtractAudio',
        'preferredcodec': 'mp3',
        'preferredquality': '192',
    }],
    'outtmpl': 'destination/%(title)s.%(ext)s',
}

with YoutubeDL(ydl_opts) as ydl:
    ydl.download([playlist_url])

#to download files from the Google colab
from google.colab import files
files.download("name of the file")
