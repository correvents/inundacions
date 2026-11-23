"""Envia un avís per Telegram i/o correu, segons els secrets que hi hagi a l'entorn.

Telegram: TELEGRAM_BOT_TOKEN i TELEGRAM_CHAT_ID.
Correu:   SMTP_USER, SMTP_PASS i AVIS_EMAIL (opcionals SMTP_HOST, per defecte smtp.gmail.com, i SMTP_PORT, 465).
Si no n'hi ha cap, només ho escriu al registre. Mai fa fallar qui el crida.

Ús des de la línia d'ordres (per provar-ho):  python3 scripts/avisa.py "Missatge de prova"
"""
import os, sys, smtplib
from email.message import EmailMessage
import requests


def envia(text, assumpte="Inundacions"):
    enviat = []
    tok, chat = os.environ.get("TELEGRAM_BOT_TOKEN"), os.environ.get("TELEGRAM_CHAT_ID")
    if tok and chat:
        try:
            r = requests.post(f"https://api.telegram.org/bot{tok}/sendMessage", timeout=30,
                              json={"chat_id": chat, "text": f"{assumpte}\n{text}", "disable_web_page_preview": True})
            r.raise_for_status()
            enviat.append("Telegram")
        except Exception as e:
            print(f"Avís per Telegram ha fallat: {type(e).__name__}", file=sys.stderr)
    usr, pwd, a = os.environ.get("SMTP_USER"), os.environ.get("SMTP_PASS"), os.environ.get("AVIS_EMAIL")
    if usr and pwd and a:
        try:
            m = EmailMessage()
            m["Subject"], m["From"], m["To"] = assumpte, usr, a
            m.set_content(text)
            with smtplib.SMTP_SSL(os.environ.get("SMTP_HOST") or "smtp.gmail.com", int(os.environ.get("SMTP_PORT") or 465), timeout=30) as s:
                s.login(usr, pwd)
                s.send_message(m)
            enviat.append("correu")
        except Exception as e:
            print(f"Avís per correu ha fallat: {type(e).__name__}", file=sys.stderr)
    print(f"Avís ({', '.join(enviat) or 'sense canal configurat'}): {assumpte} · {text}")
    return enviat


if __name__ == "__main__":
    envia(" ".join(sys.argv[1:]) or "Prova d'avís", "Inundacions · prova")
