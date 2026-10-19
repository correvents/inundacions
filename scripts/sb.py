"""Accés mínim a Supabase (PostgREST) per als scripts: esquema inundacions, clau secreta de l'entorn."""
import os, sys
import requests


def _cap(extra=None):
    clau = os.environ["SUPABASE_SECRET_KEY"]
    h = {"apikey": clau, "Authorization": f"Bearer {clau}", "Content-Profile": "inundacions", "Accept-Profile": "inundacions"}
    h.update(extra or {})
    return h


def _url(taula):
    return os.environ["SUPABASE_URL"].rstrip("/") + "/rest/v1/" + taula


def insereix(taula, files, actualitza=False):
    """Upsert per blocs. actualitza=False: els duplicats (mateixa clau primària) s'ignoren."""
    pref = "resolution=merge-duplicates" if actualitza else "resolution=ignore-duplicates"
    for i in range(0, len(files), 500):
        r = requests.post(_url(taula), json=files[i:i + 500], timeout=60, headers=_cap({"Prefer": pref + ",return=minimal"}))
        if not r.ok:
            sys.exit(f"Supabase {taula}: HTTP {r.status_code} {r.text[:300]}")
    return len(files)


def llegeix(taula, params):
    r = requests.get(_url(taula), params=params, timeout=60, headers=_cap())
    if not r.ok:
        sys.exit(f"Supabase {taula}: HTTP {r.status_code} {r.text[:300]}")
    return r.json()
